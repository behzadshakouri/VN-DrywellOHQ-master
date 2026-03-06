#!/usr/bin/env python3
"""
sip_global_equations.py

Build site-level and campaign-level global equations from the outputs of
sip_pipeline.py.

Main idea
---------
- Search a top folder recursively for:
      Results/joined__<MODE>_<TAG>Hz.csv
- Read those per-run joined tables
- Infer site name from run folder name
    Examples:
      Rosemead_S1   -> Rosemead
      Rosemead_S2   -> Rosemead
      VN_S12        -> VN
      VanNuys_S5    -> VanNuys
- Aggregate all points for each site
- Fit global equations for requested x/y combinations
- Write:
    * one combined CSV of all stacked rows
    * one Excel workbook with:
        - stacked_data_ALL
        - fits_ALL
        - settings
        - one stacked_data_<MODE> sheet per mode
        - one fits_<MODE> sheet per mode
        - one stacked_data_<MODE>_<SITE> sheet per site/mode
        - one fits_<MODE>_<SITE> sheet per site/mode
    * one plot per site/response pair
    * one overall "ALL_SITES" plot per response pair

Supported fit methods
---------------------
- linear           y = a*x + b
- origin           y = a*x
- poly             y = c0*x^d + c1*x^(d-1) + ... + cd
- logistic4        y = c + L/(1 + exp(-k*(x - x0)))

Default use case
----------------
For each site, fit:
- resistivity vs saturation
- conductivity vs saturation
- log resistivity vs log saturation
- log conductivity vs log saturation

Notes
-----
1) This script uses the joined CSV outputs from the previous pipeline.
2) "Moisture content" is interpreted as saturation S unless a theta column exists.
3) Plotting uses gnuplot (not matplotlib), so it avoids NumPy/matplotlib binary mismatch issues.
4) Plot titles/labels are plain text only; no subscript conversion.
"""

import argparse
import math
import re
import time
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================

DEFAULT_RESPONSES = [
    ("S",    "rho_real"),
    ("S",    "sig_real"),
    ("logS", "logRho"),
    ("logS", "logSig"),
]

K_MIN = 0.05
K_MAX = 50.0

# =============================================================================
# BASIC HELPERS
# =============================================================================

def fmt(x, nd=6):
    try:
        x = float(x)
        if not np.isfinite(x):
            return "nan"
        return f"{x:.{nd}g}"
    except Exception:
        return "nan"

def safe_float_series(s):
    return pd.to_numeric(s, errors="coerce").astype(float)

def sse(y, yhat):
    r = y - yhat
    return float(np.sum(r * r))

def safe_log10(x):
    x = float(x)
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return math.log10(x)

def aic_gaussian(ss_res, n, k_params):
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)

def tag_from_freq(freq: float) -> str:
    return str(freq).replace(".", "p")

def make_safe_name(s: str) -> str:
    s = str(s).strip()
    s = re.sub(r"[^\w\-\.]+", "_", s)
    return s.strip("_")

def gp_escape(s: str) -> str:
    return str(s).replace("\\", "\\\\").replace("'", "\\'")

def compact_text(s: str) -> str:
    return re.sub(r"\s+", " ", str(s)).strip()

def human_label(col: str) -> str:
    """
    Plain-text labels only. No subscript syntax, no enhanced text tricks.
    """
    mapping = {
        "S": "Saturation",
        "theta": "Moisture content",
        "rho_real": "Resistivity",
        "sig_real": "Conductivity",
        "logS": "log10(Saturation)",
        "logTheta": "log10(Moisture content)",
        "logRho": "log10(Resistivity)",
        "logSig": "log10(Conductivity)",
        "logRhoIdx": "log10(Resistivity Index)",
        "logICI_vn": "log10(ICI)",
        "logRI_vn": "log10(RI)",
        "rho_idx": "Resistivity Index",
        "sig_idx": "Conductivity Index",
    }
    return mapping.get(col, col)

def sheet_safe_name(name: str, max_len: int = 31) -> str:
    bad = r'[:\\/?*\[\]]'
    s = re.sub(bad, "_", str(name))
    return s[:max_len]

def best_fit_by_aic(fits):
    if not fits:
        return None
    good = [f for f in fits if np.isfinite(f.get("aic", np.nan))]
    if good:
        return min(good, key=lambda d: d["aic"])
    return fits[0]

def eq_with_metrics(fit):
    eq = compact_text(fit.get("equation", ""))
    parts = [eq]
    if np.isfinite(fit.get("r2", np.nan)):
        parts.append(f"R2={fmt(fit['r2'], 5)}")
    if np.isfinite(fit.get("sse", np.nan)):
        parts.append(f"SSE={fmt(fit['sse'], 5)}")
    if np.isfinite(fit.get("aic", np.nan)):
        parts.append(f"AIC={fmt(fit['aic'], 5)}")
    return " | ".join(parts)

# =============================================================================
# SITE NAME INFERENCE
# =============================================================================

def infer_site_name(run_name: str) -> str:
    s = str(run_name).strip()

    for pat in [
        r"^([A-Za-z]+(?:[A-Za-z0-9]*))(?:[_\- ]+S\d+.*)$",
        r"^([A-Za-z]+(?:[A-Za-z0-9]*))(?:[_\- ]+\d+.*)$",
        r"^([A-Za-z]+(?:[A-Za-z0-9]*))[_\- ].*$",
    ]:
        m = re.match(pat, s)
        if m:
            return m.group(1)

    return s

# =============================================================================
# FIT ROUTINES
# =============================================================================

def linear_fit(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if len(x) < 2:
        return None
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    yhat = a * x + b
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    aic = aic_gaussian(ss_res, len(x), 2)
    sign = "+" if b >= 0 else "-"
    return {
        "method": "linear",
        "n": int(len(x)),
        "a": float(a),
        "b": float(b),
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "sse": float(ss_res),
        "aic": float(aic) if np.isfinite(aic) else np.nan,
        "equation": f"y = {fmt(a)}*x {sign} {fmt(abs(b))}",
        "predict": lambda xx: a * np.asarray(xx, float) + b,
    }

def origin_fit(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if len(x) < 2:
        return None
    den = float(np.dot(x, x))
    if den <= 0:
        return None
    a = float(np.dot(x, y) / den)
    yhat = a * x
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    aic = aic_gaussian(ss_res, len(x), 1)
    return {
        "method": "origin",
        "n": int(len(x)),
        "a": float(a),
        "b": 0.0,
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "sse": float(ss_res),
        "aic": float(aic) if np.isfinite(aic) else np.nan,
        "equation": f"y = {fmt(a)}*x",
        "predict": lambda xx: a * np.asarray(xx, float),
    }

def poly_fit(x, y, deg):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if len(x) < (deg + 1):
        return None
    try:
        coeffs = np.polyfit(x, y, deg)
        yhat = np.polyval(coeffs, x)
        ss_res = sse(y, yhat)
        ss_tot = sse(y, np.full_like(y, np.mean(y)))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        aic = aic_gaussian(ss_res, len(x), deg + 1)

        terms = []
        d = len(coeffs) - 1
        for i, c in enumerate(coeffs):
            p = d - i
            cc = float(c)
            if p == 0:
                term = f"{fmt(abs(cc))}"
            elif p == 1:
                term = f"{fmt(abs(cc))}*x"
            else:
                term = f"{fmt(abs(cc))}*x^{p}"

            if not terms:
                terms.append(f"-{term}" if cc < 0 else term)
            else:
                terms.append((" - " if cc < 0 else " + ") + term)

        eq = "y = " + "".join(terms)

        return {
            "method": f"poly{deg}",
            "n": int(len(x)),
            "degree": int(deg),
            "coeffs": coeffs.tolist(),
            "r2": float(r2) if np.isfinite(r2) else np.nan,
            "sse": float(ss_res),
            "aic": float(aic) if np.isfinite(aic) else np.nan,
            "equation": eq,
            "predict": lambda xx: np.polyval(coeffs, np.asarray(xx, float)),
        }
    except Exception:
        return None

def logistic4(x, L, k, x0, c):
    z = k * (x - x0)
    z = np.clip(z, -60.0, 60.0)
    return c + L / (1.0 + np.exp(-z))

def _logistic_init_guess(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]

    y_min = float(np.min(y))
    y_max = float(np.max(y))
    L = max(1e-6, y_max - y_min)
    c = y_min

    y_mid = y_min + 0.5 * L
    j = int(np.argmin(np.abs(y - y_mid)))
    x0 = float(x[j])

    if len(x) >= 3:
        j0 = max(1, min(len(x) - 2, j))
        dy = float(y[j0 + 1] - y[j0 - 1])
        dx = float(x[j0 + 1] - x[j0 - 1])
        slope = dy / dx if dx != 0 else 0.0
        k = abs(4.0 * slope / max(L, 1e-6))
        k = float(np.clip(k, K_MIN, K_MAX))
    else:
        k = 1.0

    return L, k, x0, c

def logistic_fit_noscipy(x, y, max_seconds=8.0):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    n = len(x)

    if n < 4:
        return None

    t0 = time.time()
    Lb, kb, x0b, cb = _logistic_init_guess(x, y)

    x_min, x_max = float(np.min(x)), float(np.max(x))
    spanx = max(1e-6, x_max - x_min)

    stepL  = 0.25 * max(abs(Lb), 1e-6)
    stepk  = 0.25 * max(abs(kb), 1e-6)
    stepx0 = 0.20 * spanx
    stepc  = 0.25 * max(abs(Lb), 1e-6)

    def clamp_params(L, k, x0, c):
        L = max(1e-9, L)
        k = float(np.clip(k, K_MIN, K_MAX))
        return L, k, x0, c

    def current_sse(L, k, x0, c):
        return sse(y, logistic4(x, L, k, x0, c))

    Lb, kb, x0b, cb = clamp_params(Lb, kb, x0b, cb)
    best = current_sse(Lb, kb, x0b, cb)

    for _ in range(10):
        if max_seconds and (time.time() - t0) > max_seconds:
            break
        improved = False
        for dL in [0, -stepL, stepL]:
            for dk in [0, -stepk, stepk]:
                for dx0 in [0, -stepx0, stepx0]:
                    for dc in [0, -stepc, stepc]:
                        if dL == dk == dx0 == dc == 0:
                            continue
                        if max_seconds and (time.time() - t0) > max_seconds:
                            break
                        L2, k2, x02, c2 = clamp_params(Lb + dL, kb + dk, x0b + dx0, cb + dc)
                        s2 = current_sse(L2, k2, x02, c2)
                        if s2 + 1e-12 < best:
                            best = s2
                            Lb, kb, x0b, cb = L2, k2, x02, c2
                            improved = True
        stepL *= 0.6
        stepk *= 0.6
        stepx0 *= 0.6
        stepc *= 0.6
        if not improved and max(stepL, stepk, stepx0, stepc) < 1e-4:
            break

    yhat = logistic4(x, Lb, kb, x0b, cb)
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    aic = aic_gaussian(ss_res, n, 4)

    eq = (
        f"y = {fmt(cb)} + {fmt(Lb)}/"
        f"(1 + exp(-{fmt(kb)}*(x {'-' if x0b >= 0 else '+'} {fmt(abs(x0b))})))"
    )

    return {
        "method": "logistic4",
        "n": int(n),
        "L": float(Lb),
        "k": float(kb),
        "x0": float(x0b),
        "c": float(cb),
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "sse": float(ss_res),
        "aic": float(aic) if np.isfinite(aic) else np.nan,
        "equation": eq,
        "predict": lambda xx: logistic4(np.asarray(xx, float), Lb, kb, x0b, cb),
    }

# =============================================================================
# DATA DISCOVERY / LOADING
# =============================================================================

def find_joined_files(top: Path, mode: str, tag: str):
    pattern = f"joined__{mode}_{tag}Hz.csv"
    return sorted(top.rglob(pattern))

def add_real_alias_columns(df: pd.DataFrame):
    """
    Backfill aliases if older/newer joined files differ a bit.
    """
    out = df.copy()

    if "rho_real" not in out.columns:
        if "rho_idx" in out.columns:
            out["rho_real"] = safe_float_series(out["rho_idx"])
        elif "rho" in out.columns:
            out["rho_real"] = safe_float_series(out["rho"])

    if "sig_real" not in out.columns:
        if "sig_idx" in out.columns:
            out["sig_real"] = safe_float_series(out["sig_idx"])
        elif "sigma_imag" in out.columns:
            out["sig_real"] = safe_float_series(out["sigma_imag"])

    # Moisture-content aliases only if already present in data
    if "logTheta" not in out.columns and "theta" in out.columns:
        out["logTheta"] = out["theta"].apply(safe_log10)

    return out

def load_all_joined(top: Path, mode: str, tag: str):
    files = find_joined_files(top, mode, tag)
    rows = []

    for fp in files:
        try:
            df = pd.read_csv(fp)
        except Exception as e:
            print(f"[WARN] failed to read {fp}: {e}", flush=True)
            continue

        if df is None or df.empty:
            continue

        df = add_real_alias_columns(df)

        run_root = fp.parent.parent
        run_name = run_root.name
        site = infer_site_name(run_name)

        df["__joined_file"] = str(fp)
        df["__run_root"] = str(run_root)
        df["__run_name"] = run_name
        df["__site"] = site
        df["__mode"] = mode

        rows.append(df)

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, axis=0, ignore_index=True)

# =============================================================================
# GNUPLOT-BASED PLOTTING
# =============================================================================

def write_plot_data_csv(df_site, xcol, ycol, fits, csv_path: Path):
    x = safe_float_series(df_site[xcol]).to_numpy(float)
    y = safe_float_series(df_site[ycol]).to_numpy(float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    if len(x) == 0:
        return False

    order = np.argsort(x)
    xs = x[order]
    ys = y[order]

    xgrid = np.linspace(float(np.min(xs)), float(np.max(xs)), 300)

    n = max(len(xs), len(xgrid))
    out = pd.DataFrame(index=np.arange(n))
    out["x_data"] = np.nan
    out["y_data"] = np.nan
    out.loc[:len(xs)-1, "x_data"] = xs
    out.loc[:len(xs)-1, "y_data"] = ys

    out["x_fit"] = np.nan
    out.loc[:len(xgrid)-1, "x_fit"] = xgrid

    for fit in fits:
        col = f"y_{fit['method']}"
        out[col] = np.nan
        try:
            yfit = fit["predict"](xgrid)
            out.loc[:len(xgrid)-1, col] = yfit
        except Exception:
            pass

    out.to_csv(csv_path, index=False, na_rep="")
    return True

def write_gnuplot_script(csv_path: Path, gp_path: Path, png_path: Path,
                         title: str, xlab: str, ylab: str, fits,
                         plot_equations: bool = False,
                         eq_position: str = "right"):
    lines = [
        f"'{csv_path.name}' using (column(\"x_data\")):(column(\"y_data\")) "
        "with points pt 7 ps 1.6 title 'data'"
    ]

    for fit in fits:
        col = f"y_{fit['method']}"
        lines.append(
            f"'{csv_path.name}' using (column(\"x_fit\")):(column(\"{col}\")) "
            f"with lines lw 3 title '{gp_escape(fit['method'])}'"
        )

    plot_cmd = "plot " + ", \\\n     ".join(lines)

    eq_lines = []
    unset_lines = []
    margin_cmd = "set rmargin 10"

    if plot_equations and fits:
        if eq_position == "right":
            margin_cmd = "set rmargin 42"
            y0 = 0.94
            step = 0.08
            for i, fit in enumerate(fits[:6], start=1):
                lid = 100 + i
                txt = eq_with_metrics(fit)
                eq_lines.append(
                    f"set label {lid} '{gp_escape(txt)}' at graph 1.02,{y0:.2f} left front"
                )
                unset_lines.append(f"unset label {lid}")
                y0 -= step
        else:
            y0 = 0.95
            step = 0.07
            for i, fit in enumerate(fits[:6], start=1):
                lid = 100 + i
                txt = eq_with_metrics(fit)
                eq_lines.append(
                    f"set label {lid} '{gp_escape(txt)}' at graph 0.02,{y0:.2f} front"
                )
                unset_lines.append(f"unset label {lid}")
                y0 -= step

    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 noenhanced font 'Arial,24'
set output '{png_path.name}'

set grid
set key top right
set tics out
{margin_cmd}

set title '{gp_escape(title)}'
set xlabel '{gp_escape(xlab)}'
set ylabel '{gp_escape(ylab)}'

{chr(10).join(eq_lines)}
{plot_cmd}
{chr(10).join(unset_lines)}
"""
    gp_path.write_text(gp)

def run_gnuplot(gp_path: Path):
    subprocess.run(["gnuplot", gp_path.name], cwd=str(gp_path.parent), check=True)

def save_site_plot(df_site, xcol, ycol, fits, out_png, title,
                   plot_equations: bool = False,
                   eq_position: str = "right"):
    out_png = Path(out_png)
    csv_path = out_png.with_suffix(".plotdata.csv")
    gp_path  = out_png.with_suffix(".gp")

    ok = write_plot_data_csv(df_site, xcol, ycol, fits, csv_path)
    if not ok:
        return

    write_gnuplot_script(
        csv_path=csv_path,
        gp_path=gp_path,
        png_path=out_png,
        title=title,
        xlab=human_label(xcol),
        ylab=human_label(ycol),
        fits=fits,
        plot_equations=plot_equations,
        eq_position=eq_position,
    )
    run_gnuplot(gp_path)

# =============================================================================
# MAIN FIT DRIVER
# =============================================================================

def fit_one_group(df_group, group_name, xcol, ycol, methods, poly_deg, logistic_max_seconds):
    x = safe_float_series(df_group.get(xcol, pd.Series(dtype=float))).to_numpy(float)
    y = safe_float_series(df_group.get(ycol, pd.Series(dtype=float))).to_numpy(float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    results = []
    if len(x) < 2:
        return results

    if "linear" in methods:
        r = linear_fit(x, y)
        if r is not None:
            results.append(r)

    if "origin" in methods:
        r = origin_fit(x, y)
        if r is not None:
            results.append(r)

    if "poly" in methods and poly_deg >= 2:
        r = poly_fit(x, y, poly_deg)
        if r is not None:
            results.append(r)

    if "logistic4" in methods:
        r = logistic_fit_noscipy(x, y, max_seconds=logistic_max_seconds)
        if r is not None:
            results.append(r)

    for r in results:
        r["group"] = group_name
        r["xcol"] = xcol
        r["ycol"] = ycol

    return results

# =============================================================================
# WORKBOOK WRITING
# =============================================================================

def write_mode_site_sheets(writer, df_mode, fits_mode, mode_name: str):
    # Per-mode full sheets
    df_mode.to_excel(writer, sheet_name=sheet_safe_name(f"stacked_{mode_name}"), index=False)
    fits_mode.to_excel(writer, sheet_name=sheet_safe_name(f"fits_{mode_name}"), index=False)

    # Per-site sheets
    sites = sorted(df_mode["__site"].dropna().astype(str).unique().tolist())
    for site in sites:
        df_site = df_mode.loc[df_mode["__site"].astype(str) == site].copy()
        fits_site = fits_mode.loc[fits_mode["group"].astype(str) == site].copy() if not fits_mode.empty else pd.DataFrame()
        df_site.to_excel(
            writer,
            sheet_name=sheet_safe_name(f"data_{mode_name}_{site}"),
            index=False
        )
        if not fits_site.empty:
            fits_site.to_excel(
                writer,
                sheet_name=sheet_safe_name(f"fits_{mode_name}_{site}"),
                index=False
            )

# =============================================================================
# CLI
# =============================================================================

def parse_response_pairs(vals):
    """
    Input form:
      --pair S rho_real --pair logS logRho
    """
    pairs = []
    for v in vals:
        if len(v) != 2:
            continue
        pairs.append((v[0], v[1]))
    return pairs

def normalize_modes(mode_arg: str):
    if mode_arg == "BOTH":
        return ["NORM", "UNORM"]
    return [mode_arg]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("top", help="Top folder that contains run folders processed by sip_pipeline.py")

    ap.add_argument("--freq", type=float, default=0.01, help="Frequency tag used by previous pipeline")
    ap.add_argument("--mode", type=str, default="NORM", choices=["NORM", "UNORM", "BOTH"],
                    help="Which joined__<mode>_<tag>Hz.csv files to use")

    ap.add_argument("--pair", nargs=2, action="append", default=None,
                    help="x y pair to fit; may be repeated. Example: --pair S rho_real --pair logS logRho")

    ap.add_argument("--include-default-pairs", action="store_true",
                    help="If custom --pair values are provided, also include the default response pairs.")

    ap.add_argument("--include-vn-targets", action="store_true",
                    help="Add VN-style target pairs: logS->logRhoIdx and logS->logICI_vn")

    ap.add_argument("--include-theta", action="store_true",
                    help="If theta/logTheta exist in joined files, also add theta-based pairs automatically.")

    ap.add_argument("--methods", nargs="+", default=["linear", "poly", "logistic4"],
                    choices=["linear", "origin", "poly", "logistic4"],
                    help="Fit methods to run")

    ap.add_argument("--poly-deg", type=int, default=2, help="Polynomial degree for --methods poly")
    ap.add_argument("--logistic-max-seconds", type=float, default=20.0,
                    help="Time cap per logistic fit")

    ap.add_argument("--min-points", type=int, default=4,
                    help="Minimum points per site/response to attempt fitting")

    ap.add_argument("--site-map", type=str, default="",
                    help="Optional CSV mapping run_name,site to override inferred site names")

    ap.add_argument("--outdir", type=str, default="GLOBAL_EQUATIONS",
                    help="Output folder name created under top")

    ap.add_argument("--plot-equations", action="store_true",
                    help="Write fitted equations + metrics on global plots.")

    ap.add_argument("--eq-position", type=str, default="right", choices=["right", "inside"],
                    help="Where to place plot equations if --plot-equations is enabled.")

    ap.add_argument("--plot-best-only", action="store_true",
                    help="Plot only the best fit (by AIC) plus data, instead of all fitted curves.")

    ap.add_argument("--no-plots", action="store_true",
                    help="Skip plot generation entirely.")

    args = ap.parse_args()

    top = Path(args.top).expanduser().resolve()
    if not top.exists():
        print(f"Not found: {top}")
        raise SystemExit(2)

    outdir = top / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    tag = tag_from_freq(args.freq)
    modes = normalize_modes(args.mode)

    all_mode_frames = []
    for mode in modes:
        dfm = load_all_joined(top, mode, tag)
        if dfm.empty:
            print(f"[WARN] No joined__{mode}_{tag}Hz.csv files found under {top}", flush=True)
        else:
            all_mode_frames.append(dfm)

    if not all_mode_frames:
        print(f"No joined files found for requested mode(s) under {top}")
        raise SystemExit(1)

    df_all = pd.concat(all_mode_frames, axis=0, ignore_index=True)

    if args.site_map:
        smap = pd.read_csv(args.site_map)
        if "run_name" in smap.columns and "site" in smap.columns:
            dmap = dict(zip(smap["run_name"].astype(str), smap["site"].astype(str)))
            df_all["__site"] = df_all["__run_name"].astype(str).map(
                lambda s: dmap.get(s, infer_site_name(s))
            )

    # Build response pairs
    pairs = []
    if args.pair:
        pairs.extend(parse_response_pairs(args.pair))
        if args.include_default_pairs:
            pairs.extend(DEFAULT_RESPONSES)
    else:
        pairs.extend(DEFAULT_RESPONSES)

    if args.include_vn_targets:
        pairs.extend([
            ("logS", "logRhoIdx"),
            ("logS", "logICI_vn"),
        ])

    if args.include_theta:
        if "theta" in df_all.columns:
            pairs.extend([
                ("theta", "rho_real"),
                ("theta", "sig_real"),
            ])
        if "logTheta" in df_all.columns:
            pairs.extend([
                ("logTheta", "logRho"),
                ("logTheta", "logSig"),
            ])

    # De-duplicate pairs while preserving order
    seen = set()
    uniq_pairs = []
    for p in pairs:
        if p not in seen:
            uniq_pairs.append(p)
            seen.add(p)
    pairs = uniq_pairs

    # Save stacked data ALL
    stacked_csv = outdir / f"stacked_joined__ALLMODES_{tag}Hz.csv"
    df_all.to_csv(stacked_csv, index=False)

    fit_rows = []
    plot_jobs = []

    for mode in modes:
        dmode = df_all.loc[df_all["__mode"].astype(str) == mode].copy()
        if dmode.empty:
            continue

        groups = sorted(dmode["__site"].dropna().astype(str).unique().tolist())
        groups_all = groups + ["ALL_SITES"]

        for group_name in groups_all:
            if group_name == "ALL_SITES":
                dfg = dmode.copy()
            else:
                dfg = dmode.loc[dmode["__site"].astype(str) == group_name].copy()

            if dfg.empty:
                continue

            for xcol, ycol in pairs:
                if xcol not in dfg.columns or ycol not in dfg.columns:
                    continue

                x = safe_float_series(dfg[xcol]).to_numpy(float)
                y = safe_float_series(dfg[ycol]).to_numpy(float)
                m = np.isfinite(x) & np.isfinite(y)
                n_ok = int(np.sum(m))

                if n_ok < args.min_points:
                    continue

                fits = fit_one_group(
                    dfg, group_name, xcol, ycol,
                    methods=args.methods,
                    poly_deg=args.poly_deg,
                    logistic_max_seconds=args.logistic_max_seconds
                )

                for fr in fits:
                    row = {k: v for k, v in fr.items() if k != "predict"}
                    row["valid_points"] = n_ok
                    row["mode"] = mode
                    row["freq"] = args.freq
                    row["site"] = group_name
                    fit_rows.append(row)

                if fits and not args.no_plots:
                    fits_to_plot = [best_fit_by_aic(fits)] if args.plot_best_only else fits
                    fits_to_plot = [f for f in fits_to_plot if f is not None]
                    if fits_to_plot:
                        png = outdir / f"{make_safe_name(group_name)}__{ycol}_vs_{xcol}__{mode}_{tag}Hz.png"
                        title = f"{group_name}: global {human_label(ycol)} vs {human_label(xcol)} [{mode}, {args.freq:g} Hz]"
                        plot_jobs.append((dfg.copy(), xcol, ycol, fits_to_plot, png, title))

    for dfg, xcol, ycol, fits, png, title in plot_jobs:
        try:
            save_site_plot(
                dfg, xcol, ycol, fits, png, title,
                plot_equations=args.plot_equations,
                eq_position=args.eq_position
            )
        except Exception as e:
            print(f"[WARN] plot failed {Path(png).name}: {e}", flush=True)

    fits_df = pd.DataFrame(fit_rows)

    settings_df = pd.DataFrame({
        "key": [
            "top",
            "mode",
            "freq",
            "tag",
            "methods",
            "poly_deg",
            "logistic_max_seconds",
            "min_points",
            "response_pairs",
            "plot_equations",
            "eq_position",
            "plot_best_only",
            "no_plots",
            "include_vn_targets",
            "include_theta",
        ],
        "value": [
            str(top),
            args.mode,
            args.freq,
            tag,
            ", ".join(args.methods),
            args.poly_deg,
            args.logistic_max_seconds,
            args.min_points,
            "; ".join([f"{x}->{y}" for x, y in pairs]),
            args.plot_equations,
            args.eq_position,
            args.plot_best_only,
            args.no_plots,
            args.include_vn_targets,
            args.include_theta,
        ]
    })

    xlsx = outdir / f"global_equations__{args.mode}_{tag}Hz.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        df_all.to_excel(w, sheet_name="stacked_data_ALL", index=False)
        if not fits_df.empty:
            fits_df.to_excel(w, sheet_name="fits_ALL", index=False)
        settings_df.to_excel(w, sheet_name="settings", index=False)

        for mode in modes:
            dmode = df_all.loc[df_all["__mode"].astype(str) == mode].copy()
            fmode = fits_df.loc[fits_df["mode"].astype(str) == mode].copy() if not fits_df.empty else pd.DataFrame()
            if not dmode.empty:
                write_mode_site_sheets(w, dmode, fmode, mode)

    fits_csv = outdir / f"global_equations__{args.mode}_{tag}Hz.csv"
    fits_df.to_csv(fits_csv, index=False)

    print(f"Stacked data: {stacked_csv}", flush=True)
    print(f"Fits CSV:     {fits_csv}", flush=True)
    print(f"Fits XLSX:    {xlsx}", flush=True)
    print(f"Plots folder: {outdir}", flush=True)

if __name__ == "__main__":
    main()
