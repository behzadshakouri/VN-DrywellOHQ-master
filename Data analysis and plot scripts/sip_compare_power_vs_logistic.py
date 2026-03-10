#!/usr/bin/env python3
"""
sip_compare_power_vs_logistic.py

Compare POWER-LAW vs LOGISTIC fits using outputs previously written by sip_pipeline.py.

Reads:
    Results/joined__NORM_<tag>Hz.csv
    Results/joined__UNORM_<tag>Hz.csv

Does comparisons at 3 levels:
    1) SAMPLE   -> each run folder separately
    2) SITE     -> aggregate all runs belonging to same site
    3) ALL      -> aggregate everything

Writes a new output folder containing:
    - stacked CSVs
    - comparison CSVs
    - Excel workbooks
    - plots
    - per-mode / per-space / per-response summaries

Main goal:
    compare parameters, R2, AIC for:
        - power law
        - logistic4

Across:
    - modes: NORM, UNORM
    - spaces: REAL, LOGLOG

Priority use case:
    REAL-REAL and UNORM

Notes:
    - Plotting uses gnuplot, not matplotlib
    - Plain-text plot labels only (no subscripts)
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

K_MIN = 0.05
K_MAX = 50.0

DEFAULT_CASES = [
    # most important first
    ("UNORM", "REAL",   "S",    "rho_real"),
    ("UNORM", "REAL",   "S",    "sig_real"),
    ("NORM",  "REAL",   "S",    "rho_real"),
    ("NORM",  "REAL",   "S",    "sig_real"),
    ("UNORM", "LOGLOG", "logS", "logRho"),
    ("UNORM", "LOGLOG", "logS", "logSig"),
    ("NORM",  "LOGLOG", "logS", "logRho"),
    ("NORM",  "LOGLOG", "logS", "logSig"),
]

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

def safe_log10(x):
    x = float(x)
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return math.log10(x)

def sse(y, yhat):
    r = y - yhat
    return float(np.sum(r * r))

def aic_gaussian(ss_res, n, k_params):
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)

def tag_from_freq(freq: float) -> str:
    return str(freq).replace(".", "p")

def make_safe_name(s: str) -> str:
    s = str(s).strip()
    s = re.sub(r"[^\w\-.]+", "_", s)
    return s.strip("_")

def gp_escape(s: str) -> str:
    return str(s).replace("\\", "\\\\").replace("'", "\\'")

def compact_text(s: str) -> str:
    return re.sub(r"\s+", " ", str(s)).strip()

def human_label(col: str) -> str:
    mapping = {
        "S": "Saturation",
        "theta": "Moisture content",
        "rho_real": "Resistivity",
        "sig_real": "Conductivity",
        "rho_idx": "Resistivity Index",
        "sig_idx": "Conductivity Index",
        "logS": "log10(Saturation)",
        "logTheta": "log10(Moisture content)",
        "logRho": "log10(Resistivity)",
        "logSig": "log10(Conductivity)",
        "logRhoIdx": "log10(Resistivity Index)",
        "logICI_vn": "log10(ICI)",
        "measurement": "Measurement",
        "measurement_norm": "Measurement ID",
    }
    return mapping.get(col, col)

def sheet_safe_name(name: str, max_len: int = 31) -> str:
    bad = r'[:\\/?*\[\]]'
    s = re.sub(bad, "_", str(name))
    return s[:max_len]

# =============================================================================
# SITE / SAMPLE HELPERS
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
# DATA DISCOVERY / LOADING
# =============================================================================

def find_joined_files(top: Path, mode: str, tag: str):
    pattern = f"joined__{mode}_{tag}Hz.csv"
    return sorted(top.rglob(pattern))

def add_alias_columns(df: pd.DataFrame):
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

    if "logTheta" not in out.columns and "theta" in out.columns:
        out["logTheta"] = out["theta"].apply(safe_log10)

    return out

def load_joined_mode(top: Path, mode: str, tag: str):
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

        df = add_alias_columns(df)

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
# FITS
# =============================================================================

def r2_from_fit(y, yhat):
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return ss_res, r2

def fit_power_real(x, y):
    """
    REAL-REAL power law:
        y = a * x^b
    fitted via log-transform, evaluated in original y-space
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = x[m]
    y = y[m]
    if len(x) < 2:
        return None

    lx = np.log10(x)
    ly = np.log10(y)

    A = np.vstack([lx, np.ones_like(lx)]).T
    b, loga = np.linalg.lstsq(A, ly, rcond=None)[0]
    a = 10.0 ** float(loga)

    yhat = a * np.power(x, b)
    ss_res, r2 = r2_from_fit(y, yhat)
    aic = aic_gaussian(ss_res, len(x), 2)

    return {
        "method": "power",
        "fit_space": "REAL",
        "n": int(len(x)),
        "a": float(a),
        "b": float(b),
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "sse": float(ss_res),
        "aic": float(aic) if np.isfinite(aic) else np.nan,
        "equation": f"y = {fmt(a)} * x^{fmt(b)}",
        "predict": lambda xx: a * np.power(np.asarray(xx, float), b),
    }

def fit_power_loglog(x, y):
    """
    LOG-LOG linear form:
        y = a*x + b
    where x and y are already log-transformed columns.
    """
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
    ss_res, r2 = r2_from_fit(y, yhat)
    aic = aic_gaussian(ss_res, len(x), 2)
    sign = "+" if b >= 0 else "-"

    return {
        "method": "power",
        "fit_space": "LOGLOG",
        "n": int(len(x)),
        "a": float(a),
        "b": float(b),
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "sse": float(ss_res),
        "aic": float(aic) if np.isfinite(aic) else np.nan,
        "equation": f"y = {fmt(a)}*x {sign} {fmt(abs(b))}",
        "predict": lambda xx: a * np.asarray(xx, float) + b,
    }

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

def fit_logistic4(x, y, max_seconds=20.0):
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
    ss_res, r2 = r2_from_fit(y, yhat)
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
# COMPARISON LOGIC
# =============================================================================

def compare_two_methods(x, y, space, logistic_max_seconds):
    """
    Returns dict with both fits + winner.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    if len(x) < 4:
        return None

    if space == "REAL":
        power = fit_power_real(x, y)
    else:
        power = fit_power_loglog(x, y)

    logistic = fit_logistic4(x, y, max_seconds=logistic_max_seconds)

    if power is None or logistic is None:
        return None

    if np.isfinite(power["aic"]) and np.isfinite(logistic["aic"]):
        winner = "power" if power["aic"] < logistic["aic"] else "logistic4"
    elif np.isfinite(power["r2"]) and np.isfinite(logistic["r2"]):
        winner = "power" if power["r2"] > logistic["r2"] else "logistic4"
    else:
        winner = "unknown"

    return {
        "power": power,
        "logistic4": logistic,
        "winner": winner,
    }

def add_compare_row(rows, level, group_name, mode, space, xcol, ycol,
                    n_valid, folder_sources, run_names, fitcmp):
    p = fitcmp["power"]
    l = fitcmp["logistic4"]

    rows.append({
        "level": level,
        "group": group_name,
        "mode": mode,
        "space": space,
        "xcol": xcol,
        "ycol": ycol,
        "n_valid": n_valid,
        "winner": fitcmp["winner"],

        "power_equation": p.get("equation", ""),
        "power_a": p.get("a", np.nan),
        "power_b": p.get("b", np.nan),
        "power_r2": p.get("r2", np.nan),
        "power_aic": p.get("aic", np.nan),
        "power_sse": p.get("sse", np.nan),

        "logistic_equation": l.get("equation", ""),
        "logistic_L": l.get("L", np.nan),
        "logistic_k": l.get("k", np.nan),
        "logistic_x0": l.get("x0", np.nan),
        "logistic_c": l.get("c", np.nan),
        "logistic_r2": l.get("r2", np.nan),
        "logistic_aic": l.get("aic", np.nan),
        "logistic_sse": l.get("sse", np.nan),

        "delta_r2_logistic_minus_power": (
            l.get("r2", np.nan) - p.get("r2", np.nan)
            if np.isfinite(l.get("r2", np.nan)) and np.isfinite(p.get("r2", np.nan))
            else np.nan
        ),
        "delta_aic_logistic_minus_power": (
            l.get("aic", np.nan) - p.get("aic", np.nan)
            if np.isfinite(l.get("aic", np.nan)) and np.isfinite(p.get("aic", np.nan))
            else np.nan
        ),

        "folder_sources": folder_sources,
        "run_names": run_names,
    })

# =============================================================================
# PLOTTING
# =============================================================================

def write_plot_data_csv(df_used, xcol, ycol, fitcmp, csv_path: Path):
    x = safe_float_series(df_used[xcol]).to_numpy(float)
    y = safe_float_series(df_used[ycol]).to_numpy(float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    if len(x) == 0:
        return False

    order = np.argsort(x)
    xs = x[order]
    ys = y[order]

    xgrid = np.linspace(float(np.min(xs)), float(np.max(xs)), 300)

    power = fitcmp["power"]
    logistic = fitcmp["logistic4"]

    n = max(len(xs), len(xgrid))
    out = pd.DataFrame(index=np.arange(n))
    out["x_data"] = np.nan
    out["y_data"] = np.nan
    out.loc[:len(xs)-1, "x_data"] = xs
    out.loc[:len(xs)-1, "y_data"] = ys

    out["x_fit"] = np.nan
    out.loc[:len(xgrid)-1, "x_fit"] = xgrid

    out["y_power"] = np.nan
    out["y_logistic4"] = np.nan
    out.loc[:len(xgrid)-1, "y_power"] = power["predict"](xgrid)
    out.loc[:len(xgrid)-1, "y_logistic4"] = logistic["predict"](xgrid)

    out.to_csv(csv_path, index=False, na_rep="")
    return True

def write_gnuplot_script(csv_path: Path, gp_path: Path, png_path: Path,
                         title: str, xlab: str, ylab: str,
                         fitcmp, plot_equations: bool = False,
                         eq_position: str = "right"):
    eq_lines = []
    unset_lines = []
    rmargin = "set rmargin 10"

    ptxt = (
        f"POWER: {compact_text(fitcmp['power']['equation'])} | "
        f"R2={fmt(fitcmp['power']['r2'],5)} | AIC={fmt(fitcmp['power']['aic'],5)}"
    )
    ltxt = (
        f"LOGISTIC: {compact_text(fitcmp['logistic4']['equation'])} | "
        f"R2={fmt(fitcmp['logistic4']['r2'],5)} | AIC={fmt(fitcmp['logistic4']['aic'],5)}"
    )
    wtxt = f"Winner: {fitcmp['winner']}"

    if plot_equations:
        if eq_position == "right":
            rmargin = "set rmargin 44"
            eq_lines.extend([
                f"set label 101 '{gp_escape(ptxt)}' at graph 1.02,0.94 left front",
                f"set label 102 '{gp_escape(ltxt)}' at graph 1.02,0.84 left front",
                f"set label 103 '{gp_escape(wtxt)}' at graph 1.02,0.74 left front",
            ])
        else:
            eq_lines.extend([
                f"set label 101 '{gp_escape(ptxt)}' at graph 0.02,0.95 front",
                f"set label 102 '{gp_escape(ltxt)}' at graph 0.02,0.87 front",
                f"set label 103 '{gp_escape(wtxt)}' at graph 0.02,0.79 front",
            ])
        unset_lines.extend(["unset label 101", "unset label 102", "unset label 103"])

    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 noenhanced font 'Arial,24'
set output '{png_path.name}'

set grid
set key top right
set tics out
{rmargin}

set title '{gp_escape(title)}'
set xlabel '{gp_escape(xlab)}'
set ylabel '{gp_escape(ylab)}'

{chr(10).join(eq_lines)}

plot '{csv_path.name}' using (column("x_data")):(column("y_data")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_path.name}' using (column("x_fit")):(column("y_power")) with lines lw 3 title 'power', \\
     '{csv_path.name}' using (column("x_fit")):(column("y_logistic4")) with lines lw 3 title 'logistic4'

{chr(10).join(unset_lines)}
"""
    gp_path.write_text(gp)

def run_gnuplot(gp_path: Path):
    subprocess.run(["gnuplot", gp_path.name], cwd=str(gp_path.parent), check=True)

def save_compare_plot(df_used, xcol, ycol, fitcmp, out_png, title,
                      plot_equations=False, eq_position="right"):
    out_png = Path(out_png)
    csv_path = out_png.with_suffix(".plotdata.csv")
    gp_path = out_png.with_suffix(".gp")

    ok = write_plot_data_csv(df_used, xcol, ycol, fitcmp, csv_path)
    if not ok:
        return

    write_gnuplot_script(
        csv_path=csv_path,
        gp_path=gp_path,
        png_path=out_png,
        title=title,
        xlab=human_label(xcol),
        ylab=human_label(ycol),
        fitcmp=fitcmp,
        plot_equations=plot_equations,
        eq_position=eq_position,
    )
    run_gnuplot(gp_path)

# =============================================================================
# CASE HANDLING
# =============================================================================

def normalize_modes(mode_arg: str):
    if mode_arg == "BOTH":
        return ["NORM", "UNORM"]
    return [mode_arg]

def parse_cases(custom_cases):
    """
    Input:
      --case UNORM REAL S rho_real
    """
    out = []
    for item in custom_cases or []:
        if len(item) != 4:
            continue
        out.append((item[0], item[1], item[2], item[3]))
    return out

# =============================================================================
# MAIN
# =============================================================================

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("top", help="Top folder containing previous sip_pipeline outputs")
    ap.add_argument("--freq", type=float, default=0.01, help="Frequency tag")
    ap.add_argument("--mode", type=str, default="BOTH", choices=["NORM", "UNORM", "BOTH"],
                    help="Which mode(s) to process")

    ap.add_argument("--case", nargs=4, action="append", default=None,
                    help="Custom case: MODE SPACE XCOL YCOL ; repeatable")

    ap.add_argument("--include-default-cases", action="store_true",
                    help="If custom --case values are provided, also include default cases")

    ap.add_argument("--include-vn-targets", action="store_true",
                    help="Also compare logS->logRhoIdx and logS->logICI_vn in LOGLOG")

    ap.add_argument("--include-theta", action="store_true",
                    help="Also compare theta/logTheta if present")

    ap.add_argument("--logistic-max-seconds", type=float, default=20.0,
                    help="Time cap per logistic fit")

    ap.add_argument("--min-points", type=int, default=4,
                    help="Minimum points required")

    ap.add_argument("--plot-equations", action="store_true",
                    help="Write equations and R2/AIC on plots")

    ap.add_argument("--eq-position", type=str, default="right", choices=["right", "inside"],
                    help="Equation placement on plots")

    ap.add_argument("--no-plots", action="store_true",
                    help="Skip plot generation")

    ap.add_argument("--site-map", type=str, default="",
                    help="Optional CSV mapping run_name,site to override inferred site names")

    ap.add_argument("--outdir", type=str, default="COMPARE_POWER_VS_LOGISTIC",
                    help="Output folder under top")

    args = ap.parse_args()

    top = Path(args.top).expanduser().resolve()
    if not top.exists():
        print(f"Not found: {top}")
        raise SystemExit(2)

    outdir = top / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    tag = tag_from_freq(args.freq)
    modes = normalize_modes(args.mode)

    # load data
    frames = []
    for mode in modes:
        dfm = load_joined_mode(top, mode, tag)
        if dfm.empty:
            print(f"[WARN] no joined__{mode}_{tag}Hz.csv files found", flush=True)
        else:
            frames.append(dfm)

    if not frames:
        print("No joined files found.")
        raise SystemExit(1)

    df_all = pd.concat(frames, axis=0, ignore_index=True)

    if args.site_map:
        smap = pd.read_csv(args.site_map)
        if "run_name" in smap.columns and "site" in smap.columns:
            dmap = dict(zip(smap["run_name"].astype(str), smap["site"].astype(str)))
            df_all["__site"] = df_all["__run_name"].astype(str).map(lambda s: dmap.get(s, infer_site_name(s)))

    # cases
    cases = []
    user_cases = parse_cases(args.case)
    if user_cases:
        cases.extend(user_cases)
        if args.include_default_cases:
            cases.extend(DEFAULT_CASES)
    else:
        cases.extend(DEFAULT_CASES)

    if args.include_vn_targets:
        for m in modes:
            cases.extend([
                (m, "LOGLOG", "logS", "logRhoIdx"),
                (m, "LOGLOG", "logS", "logICI_vn"),
            ])

    if args.include_theta:
        for m in modes:
            if "theta" in df_all.columns:
                cases.extend([
                    (m, "REAL", "theta", "rho_real"),
                    (m, "REAL", "theta", "sig_real"),
                ])
            if "logTheta" in df_all.columns:
                cases.extend([
                    (m, "LOGLOG", "logTheta", "logRho"),
                    (m, "LOGLOG", "logTheta", "logSig"),
                ])

    # deduplicate cases
    seen = set()
    uniq = []
    for c in cases:
        if c not in seen:
            uniq.append(c)
            seen.add(c)
    cases = uniq

    # keep only requested modes
    cases = [c for c in cases if c[0] in modes]

    # write stacked all data
    stacked_csv = outdir / f"stacked_compare_input__{args.mode}_{tag}Hz.csv"
    df_all.to_csv(stacked_csv, index=False)

    compare_rows = []
    plot_jobs = []
    used_rows_all = []

    levels = ["SAMPLE", "SITE", "ALL"]

    for mode, space, xcol, ycol in cases:
        dmode = df_all.loc[df_all["__mode"].astype(str) == mode].copy()
        if dmode.empty:
            continue
        if xcol not in dmode.columns or ycol not in dmode.columns:
            continue

        for level in levels:
            if level == "SAMPLE":
                groups = sorted(dmode["__run_name"].astype(str).unique().tolist())
                group_key = "__run_name"
            elif level == "SITE":
                groups = sorted(dmode["__site"].astype(str).unique().tolist())
                group_key = "__site"
            else:
                groups = ["ALL_DATA"]
                group_key = None

            for g in groups:
                if level == "ALL":
                    dfg = dmode.copy()
                else:
                    dfg = dmode.loc[dmode[group_key].astype(str) == g].copy()

                if dfg.empty:
                    continue

                x = safe_float_series(dfg[xcol]).to_numpy(float)
                y = safe_float_series(dfg[ycol]).to_numpy(float)
                m = np.isfinite(x) & np.isfinite(y)
                n_valid = int(np.sum(m))
                if n_valid < args.min_points:
                    continue

                fitcmp = compare_two_methods(
                    x=x[m],
                    y=y[m],
                    space=space,
                    logistic_max_seconds=args.logistic_max_seconds,
                )
                if fitcmp is None:
                    continue

                folder_sources = "; ".join(sorted(dfg["__run_root"].astype(str).unique().tolist()))
                run_names = "; ".join(sorted(dfg["__run_name"].astype(str).unique().tolist()))

                add_compare_row(
                    rows=compare_rows,
                    level=level,
                    group_name=g,
                    mode=mode,
                    space=space,
                    xcol=xcol,
                    ycol=ycol,
                    n_valid=n_valid,
                    folder_sources=folder_sources,
                    run_names=run_names,
                    fitcmp=fitcmp,
                )

                dfg_used = dfg.loc[m].copy()
                dfg_used["__compare_level"] = level
                dfg_used["__compare_group"] = g
                dfg_used["__space"] = space
                dfg_used["__xcol"] = xcol
                dfg_used["__ycol"] = ycol
                dfg_used["__winner"] = fitcmp["winner"]
                used_rows_all.append(dfg_used)

                if not args.no_plots:
                    png = outdir / (
                        f"{make_safe_name(level)}__{make_safe_name(g)}__"
                        f"{ycol}_vs_{xcol}__{mode}_{space}_{tag}Hz.png"
                    )
                    title = (
                        f"{level} {g}: {human_label(ycol)} vs {human_label(xcol)} "
                        f"[{mode}, {space}, {args.freq:g} Hz]"
                    )
                    plot_jobs.append((dfg_used.copy(), xcol, ycol, fitcmp, png, title))

    compare_df = pd.DataFrame(compare_rows)
    used_df = pd.concat(used_rows_all, axis=0, ignore_index=True) if used_rows_all else pd.DataFrame()

    for dfg_used, xcol, ycol, fitcmp, png, title in plot_jobs:
        try:
            save_compare_plot(
                dfg_used, xcol, ycol, fitcmp, png, title,
                plot_equations=args.plot_equations,
                eq_position=args.eq_position,
            )
        except Exception as e:
            print(f"[WARN] plot failed {png.name}: {e}", flush=True)

    # outputs
    compare_csv = outdir / f"compare_power_vs_logistic__{args.mode}_{tag}Hz.csv"
    compare_df.to_csv(compare_csv, index=False)

    used_csv = outdir / f"compare_power_vs_logistic_used_data__{args.mode}_{tag}Hz.csv"
    used_df.to_csv(used_csv, index=False)

    settings_df = pd.DataFrame({
        "key": [
            "top",
            "mode",
            "freq",
            "tag",
            "logistic_max_seconds",
            "min_points",
            "plot_equations",
            "eq_position",
            "no_plots",
            "cases",
        ],
        "value": [
            str(top),
            args.mode,
            args.freq,
            tag,
            args.logistic_max_seconds,
            args.min_points,
            args.plot_equations,
            args.eq_position,
            args.no_plots,
            "; ".join([f"{m}|{sp}|{x}|{y}" for (m, sp, x, y) in cases]),
        ]
    })

    xlsx = outdir / f"compare_power_vs_logistic__{args.mode}_{tag}Hz.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        if not used_df.empty:
            used_df.to_excel(w, sheet_name="used_data_ALL", index=False)
        if not compare_df.empty:
            compare_df.to_excel(w, sheet_name="compare_ALL", index=False)
        settings_df.to_excel(w, sheet_name="settings", index=False)

        # per mode
        for mode in modes:
            dmode_used = used_df.loc[used_df["__mode"].astype(str) == mode].copy() if not used_df.empty else pd.DataFrame()
            dmode_cmp = compare_df.loc[compare_df["mode"].astype(str) == mode].copy() if not compare_df.empty else pd.DataFrame()

            if not dmode_used.empty:
                dmode_used.to_excel(w, sheet_name=sheet_safe_name(f"used_{mode}"), index=False)
            if not dmode_cmp.empty:
                dmode_cmp.to_excel(w, sheet_name=sheet_safe_name(f"compare_{mode}"), index=False)

        # per level
        for level in ["SAMPLE", "SITE", "ALL"]:
            dlev = compare_df.loc[compare_df["level"].astype(str) == level].copy() if not compare_df.empty else pd.DataFrame()
            if not dlev.empty:
                dlev.to_excel(w, sheet_name=sheet_safe_name(f"cmp_{level}"), index=False)

    # per-mode xlsx
    mode_xlsx_paths = []
    for mode in modes:
        dmode_used = used_df.loc[used_df["__mode"].astype(str) == mode].copy() if not used_df.empty else pd.DataFrame()
        dmode_cmp = compare_df.loc[compare_df["mode"].astype(str) == mode].copy() if not compare_df.empty else pd.DataFrame()

        mode_xlsx = outdir / f"compare_power_vs_logistic__{mode}_{tag}Hz.xlsx"
        with pd.ExcelWriter(mode_xlsx, engine="openpyxl") as w:
            if not dmode_used.empty:
                dmode_used.to_excel(w, sheet_name="used_data", index=False)
            if not dmode_cmp.empty:
                dmode_cmp.to_excel(w, sheet_name="compare", index=False)
            settings_df.to_excel(w, sheet_name="settings", index=False)
        mode_xlsx_paths.append(mode_xlsx)

    print(f"Used data CSV: {used_csv}", flush=True)
    print(f"Compare CSV:   {compare_csv}", flush=True)
    print(f"Compare XLSX:  {xlsx}", flush=True)
    for p in mode_xlsx_paths:
        print(f"Mode XLSX:     {p}", flush=True)
    print(f"Plots folder:  {outdir}", flush=True)

if __name__ == "__main__":
    main()
