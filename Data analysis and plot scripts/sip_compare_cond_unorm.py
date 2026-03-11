#!/usr/bin/env python3
"""
sip_compare_cond_unorm.py

Focused comparison script for CONDUCTIVITY with UNORM only.

What it does
------------
Reads:
    Results/joined__UNORM_<tag>Hz.csv

Compares POWER-LAW vs LOGISTIC4 for conductivity only:
    x = S
    y = sig_real

at three aggregation levels:
    1) SAMPLE
    2) SITE
    3) ALL

Important additions
-------------------
1) Writes, for BOTH methods:
   - residual std (sigma_hat)
   - Gaussian log-likelihood (lnL)
   - AIC
   so model choice can be made directly by smaller AIC.

2) Also writes logistic amplitude parameter separately:
   - logistic_amp_L
   This avoids confusion between:
   - likelihood L / lnL
   - logistic parameter L

3) Adds switch for logarithm base used by power-law transform:
   --log-base 10
   --log-base e

4) Normalization is done on the RELEVANT FIT RANGE ONLY:
   - SAMPLE fits use sample-specific normalization
   - SITE fits use site-specific normalization
   - ALL fits use all-data normalization
   So no sample-by-sample normalization leaks into SITE or ALL.

5) Allows monotonic decreasing logistic behavior:
   as S increases, the fitted logistic can decrease.
   This is done by allowing signed logistic amplitude.

Outputs
-------
Creates an output folder with:
    - stacked input CSV
    - used data CSV
    - comparison CSV
    - Excel workbook
    - gnuplot plots

Notes
-----
- Plotting uses gnuplot
- Plain-text labels only
- Default comparison is in REAL space:
      sig_real vs S
- Optional log-space comparison:
      log(sig_real) vs log(S)
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

K_ABS_MIN = 0.05
K_ABS_MAX = 50.0

# default cases: focused on conductivity only
DEFAULT_CASES = [
    ("REAL", "S", "sig_real"),
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

def safe_ln(x):
    x = float(x)
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return math.log(x)

def apply_log_scalar(x, base: str):
    if base == "e":
        return safe_ln(x)
    return safe_log10(x)

def apply_log_array(x, base: str):
    x = np.asarray(x, float)
    out = np.full_like(x, np.nan, dtype=float)
    m = np.isfinite(x) & (x > 0)
    if base == "e":
        out[m] = np.log(x[m])
    else:
        out[m] = np.log10(x[m])
    return out

def inverse_log_array(x, base: str):
    x = np.asarray(x, float)
    if base == "e":
        return np.exp(x)
    return np.power(10.0, x)

def sse(y, yhat):
    r = y - yhat
    return float(np.sum(r * r))

def r2_from_fit(y, yhat):
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return ss_res, r2

def sigma_hat_from_sse(ss_res, n):
    if n <= 0 or not np.isfinite(ss_res) or ss_res < 0:
        return np.nan
    return float(np.sqrt(ss_res / n))

def gaussian_loglik(ss_res, n):
    """
    Gaussian MLE log-likelihood with sigma^2 = SSE/n.
    """
    if n <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    sigma2 = ss_res / n
    return float(-0.5 * n * (math.log(2.0 * math.pi * sigma2) + 1.0))

def aic_gaussian(ss_res, n, k_params):
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    # equivalent to -2 lnL + 2k for Gaussian MLE
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

def human_label(col: str):
    mapping = {
        "S": "Saturation",
        "sig_real": "Conductivity",
        "logS": "log(Saturation)",
        "logSig": "log(Conductivity)",
    }
    return mapping.get(col, col)

def sheet_safe_name(name: str, max_len: int = 31) -> str:
    bad = r'[:\\/?*\[\]]'
    s = re.sub(bad, "_", str(name))
    return s[:max_len]

# =============================================================================
# SITE / FILE HELPERS
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

def find_joined_files(top: Path, tag: str):
    pattern = f"joined__UNORM_{tag}Hz.csv"
    return sorted(top.rglob(pattern))

def add_alias_columns(df: pd.DataFrame):
    out = df.copy()

    if "sig_real" not in out.columns:
        if "sig_idx" in out.columns:
            out["sig_real"] = safe_float_series(out["sig_idx"])
        elif "sigma_imag" in out.columns:
            out["sig_real"] = safe_float_series(out["sigma_imag"])

    if "logS" not in out.columns and "S" in out.columns:
        out["logS"] = out["S"].apply(safe_log10)

    if "logSig" not in out.columns and "sig_real" in out.columns:
        out["logSig"] = out["sig_real"].apply(safe_log10)

    return out

def load_unorm_joined(top: Path, tag: str):
    files = find_joined_files(top, tag)
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
        df["__mode"] = "UNORM"

        rows.append(df)

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, axis=0, ignore_index=True)

# =============================================================================
# NORMALIZATION
# =============================================================================

def fit_normalizer(v, method: str):
    v = np.asarray(v, float)
    m = np.isfinite(v)
    vv = v[m]

    if len(vv) == 0 or method == "none":
        return {"method": "none", "center": 0.0, "scale": 1.0, "vmin": np.nan, "vmax": np.nan}

    if method == "zscore":
        mu = float(np.mean(vv))
        sd = float(np.std(vv))
        if not np.isfinite(sd) or sd <= 0:
            sd = 1.0
        return {"method": "zscore", "center": mu, "scale": sd, "vmin": np.nan, "vmax": np.nan}

    if method == "minmax":
        vmin = float(np.min(vv))
        vmax = float(np.max(vv))
        rng = vmax - vmin
        if not np.isfinite(rng) or rng <= 0:
            rng = 1.0
        return {"method": "minmax", "center": vmin, "scale": rng, "vmin": vmin, "vmax": vmax}

    return {"method": "none", "center": 0.0, "scale": 1.0, "vmin": np.nan, "vmax": np.nan}

def apply_normalizer(v, norm):
    v = np.asarray(v, float)
    m = np.isfinite(v)
    out = np.full_like(v, np.nan, dtype=float)

    if norm["method"] == "none":
        out[m] = v[m]
        return out

    out[m] = (v[m] - norm["center"]) / norm["scale"]
    return out

def inverse_normalizer(v, norm):
    v = np.asarray(v, float)

    if norm["method"] == "none":
        return v

    return v * norm["scale"] + norm["center"]

# =============================================================================
# FITS
# =============================================================================

def fit_power_real(x, y, log_base="10"):
    """
    REAL-space power law:
        y = a * x^b
    fit by log-transform, evaluate in REAL y-space.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = x[m]
    y = y[m]
    if len(x) < 2:
        return None

    lx = apply_log_array(x, log_base)
    ly = apply_log_array(y, log_base)

    A = np.vstack([lx, np.ones_like(lx)]).T
    b, loga = np.linalg.lstsq(A, ly, rcond=None)[0]

    if log_base == "e":
        a = float(math.exp(loga))
    else:
        a = float(10.0 ** loga)

    yhat = a * np.power(x, b)
    ss_res, r2 = r2_from_fit(y, yhat)
    sig = sigma_hat_from_sse(ss_res, len(x))
    lnL = gaussian_loglik(ss_res, len(x))
    aic = aic_gaussian(ss_res, len(x), 2)

    return {
        "method": "power",
        "space": "REAL",
        "n": int(len(x)),
        "a": float(a),
        "b": float(b),
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "sse": float(ss_res),
        "sigma_hat": float(sig) if np.isfinite(sig) else np.nan,
        "lnL": float(lnL) if np.isfinite(lnL) else np.nan,
        "aic": float(aic) if np.isfinite(aic) else np.nan,
        "equation": f"y = {fmt(a)} * x^{fmt(b)}",
        "predict": lambda xx: a * np.power(np.asarray(xx, float), b),
    }

def fit_power_logspace(x, y):
    """
    LOG-space linear form:
        y = a*x + b
    where x and y are already transformed logs.
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
    sig = sigma_hat_from_sse(ss_res, len(x))
    lnL = gaussian_loglik(ss_res, len(x))
    aic = aic_gaussian(ss_res, len(x), 2)
    sign = "+" if b >= 0 else "-"

    return {
        "method": "power",
        "space": "LOG",
        "n": int(len(x)),
        "a": float(a),
        "b": float(b),
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "sse": float(ss_res),
        "sigma_hat": float(sig) if np.isfinite(sig) else np.nan,
        "lnL": float(lnL) if np.isfinite(lnL) else np.nan,
        "aic": float(aic) if np.isfinite(aic) else np.nan,
        "equation": f"y = {fmt(a)}*x {sign} {fmt(abs(b))}",
        "predict": lambda xx: a * np.asarray(xx, float) + b,
    }

def logistic4_signed(x, L, k, x0, c):
    z = k * (x - x0)
    z = np.clip(z, -60.0, 60.0)
    return c + L / (1.0 + np.exp(-z))

def _logistic_init_guess_signed(x, y):
    """
    Allows decreasing curves by allowing signed L.
    For increasing:
        c ~ low plateau, L > 0
    For decreasing:
        c ~ high plateau, L < 0
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)

    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]

    y_min = float(np.min(y))
    y_max = float(np.max(y))
    yr = max(1e-6, y_max - y_min)

    trend = float(y[-1] - y[0])

    if trend >= 0:
        c = y_min
        L = yr
        y_mid = c + 0.5 * L
    else:
        c = y_max
        L = -yr
        y_mid = c + 0.5 * L

    j = int(np.argmin(np.abs(y - y_mid)))
    x0 = float(x[j])

    if len(x) >= 3:
        j0 = max(1, min(len(x) - 2, j))
        dy = float(y[j0 + 1] - y[j0 - 1])
        dx = float(x[j0 + 1] - x[j0 - 1])
        slope = dy / dx if dx != 0 else 0.0
        k = abs(4.0 * slope / max(abs(L), 1e-6))
        k = float(np.clip(k, K_ABS_MIN, K_ABS_MAX))
    else:
        k = 1.0

    return L, k, x0, c

def fit_logistic4_signed(x, y, max_seconds=20.0):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    n = len(x)
    if n < 4:
        return None

    t0 = time.time()
    Lb, kb, x0b, cb = _logistic_init_guess_signed(x, y)

    x_min, x_max = float(np.min(x)), float(np.max(x))
    spanx = max(1e-6, x_max - x_min)

    stepL  = 0.25 * max(abs(Lb), 1e-6)
    stepk  = 0.25 * max(abs(kb), 1e-6)
    stepx0 = 0.20 * spanx
    stepc  = 0.25 * max(abs(Lb), 1e-6)

    def clamp_params(L, k, x0, c):
        # signed L allowed; k kept positive
        if abs(L) < 1e-9:
            L = -1e-9 if L < 0 else 1e-9
        k = float(np.clip(abs(k), K_ABS_MIN, K_ABS_MAX))
        return float(L), float(k), float(x0), float(c)

    def current_sse(L, k, x0, c):
        return sse(y, logistic4_signed(x, L, k, x0, c))

    Lb, kb, x0b, cb = clamp_params(Lb, kb, x0b, cb)
    best = current_sse(Lb, kb, x0b, cb)

    for _ in range(12):
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

    yhat = logistic4_signed(x, Lb, kb, x0b, cb)
    ss_res, r2 = r2_from_fit(y, yhat)
    sig = sigma_hat_from_sse(ss_res, n)
    lnL = gaussian_loglik(ss_res, n)
    aic = aic_gaussian(ss_res, n, 4)

    eq = (
        f"y = {fmt(cb)} + {fmt(Lb)}/"
        f"(1 + exp(-{fmt(kb)}*(x {'-' if x0b >= 0 else '+'} {fmt(abs(x0b))})))"
    )

    return {
        "method": "logistic4",
        "n": int(n),
        "amp_L": float(Lb),   # signed amplitude
        "k": float(kb),
        "x0": float(x0b),
        "c": float(cb),
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "sse": float(ss_res),
        "sigma_hat": float(sig) if np.isfinite(sig) else np.nan,
        "lnL": float(lnL) if np.isfinite(lnL) else np.nan,
        "aic": float(aic) if np.isfinite(aic) else np.nan,
        "equation": eq,
        "predict": lambda xx: logistic4_signed(np.asarray(xx, float), Lb, kb, x0b, cb),
    }

# =============================================================================
# COMPARISON
# =============================================================================

def choose_winner(power, logistic):
    if power is None or logistic is None:
        return "unknown"

    paic = power.get("aic", np.nan)
    laic = logistic.get("aic", np.nan)

    if np.isfinite(paic) and np.isfinite(laic):
        if paic < laic:
            return "power"
        elif laic < paic:
            return "logistic4"

    pr2 = power.get("r2", np.nan)
    lr2 = logistic.get("r2", np.nan)
    if np.isfinite(pr2) and np.isfinite(lr2):
        if pr2 > lr2:
            return "power"
        elif lr2 > pr2:
            return "logistic4"

    return "tie_or_unknown"

def compare_two_methods(x, y, space, log_base, logistic_max_seconds):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    if len(x) < 4:
        return None

    if space == "REAL":
        power = fit_power_real(x, y, log_base=log_base)
    else:
        power = fit_power_logspace(x, y)

    logistic = fit_logistic4_signed(x, y, max_seconds=logistic_max_seconds)

    if power is None or logistic is None:
        return None

    winner = choose_winner(power, logistic)

    return {
        "power": power,
        "logistic4": logistic,
        "winner": winner,
    }

# =============================================================================
# PLOTTING
# =============================================================================

def write_plot_data_csv(df_used, xcol_fit, ycol_fit, fitcmp, csv_path: Path):
    x = safe_float_series(df_used[xcol_fit]).to_numpy(float)
    y = safe_float_series(df_used[ycol_fit]).to_numpy(float)
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
        f"sigma={fmt(fitcmp['power']['sigma_hat'],5)} | "
        f"lnL={fmt(fitcmp['power']['lnL'],5)} | "
        f"AIC={fmt(fitcmp['power']['aic'],5)}"
    )
    ltxt = (
        f"LOGISTIC: {compact_text(fitcmp['logistic4']['equation'])} | "
        f"amp_L={fmt(fitcmp['logistic4']['amp_L'],5)} | "
        f"sigma={fmt(fitcmp['logistic4']['sigma_hat'],5)} | "
        f"lnL={fmt(fitcmp['logistic4']['lnL'],5)} | "
        f"AIC={fmt(fitcmp['logistic4']['aic'],5)}"
    )
    wtxt = f"Winner: {fitcmp['winner']}"

    if plot_equations:
        if eq_position == "right":
            rmargin = "set rmargin 48"
            eq_lines.extend([
                f"set label 101 '{gp_escape(ptxt)}' at graph 1.02,0.94 left front",
                f"set label 102 '{gp_escape(ltxt)}' at graph 1.02,0.82 left front",
                f"set label 103 '{gp_escape(wtxt)}' at graph 1.02,0.70 left front",
            ])
        else:
            eq_lines.extend([
                f"set label 101 '{gp_escape(ptxt)}' at graph 0.02,0.95 front",
                f"set label 102 '{gp_escape(ltxt)}' at graph 0.02,0.86 front",
                f"set label 103 '{gp_escape(wtxt)}' at graph 0.02,0.77 front",
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

def save_compare_plot(df_used, xcol_fit, ycol_fit, fitcmp, out_png, title,
                      plot_equations=False, eq_position="right"):
    out_png = Path(out_png)
    csv_path = out_png.with_suffix(".plotdata.csv")
    gp_path = out_png.with_suffix(".gp")

    ok = write_plot_data_csv(df_used, xcol_fit, ycol_fit, fitcmp, csv_path)
    if not ok:
        return

    write_gnuplot_script(
        csv_path=csv_path,
        gp_path=gp_path,
        png_path=out_png,
        title=title,
        xlab=human_label(xcol_fit),
        ylab=human_label(ycol_fit),
        fitcmp=fitcmp,
        plot_equations=plot_equations,
        eq_position=eq_position,
    )
    run_gnuplot(gp_path)

# =============================================================================
# GROUP PREP
# =============================================================================

def prepare_case_dataframe(dfg, space: str, xcol_raw: str, ycol_raw: str,
                           log_base: str, normalize: str):
    """
    Build the actual x/y used for fitting for ONE group only.
    Normalization is therefore group-specific by construction.
    """
    out = dfg.copy()

    if xcol_raw not in out.columns or ycol_raw not in out.columns:
        return None, None

    x_raw = safe_float_series(out[xcol_raw]).to_numpy(float)
    y_raw = safe_float_series(out[ycol_raw]).to_numpy(float)

    if space == "LOG":
        x_fit = apply_log_array(x_raw, log_base)
        y_fit = apply_log_array(y_raw, log_base)
        x_fit_name = "x_fit_log"
        y_fit_name = "y_fit_log"
    else:
        x_fit = x_raw.copy()
        y_fit = y_raw.copy()
        x_fit_name = "x_fit_real"
        y_fit_name = "y_fit_real"

    # mask BEFORE normalization
    m = np.isfinite(x_fit) & np.isfinite(y_fit)
    if np.sum(m) == 0:
        return None, None

    x_norm = fit_normalizer(x_fit[m], normalize)
    y_norm = fit_normalizer(y_fit[m], normalize)

    x_used = np.full_like(x_fit, np.nan, dtype=float)
    y_used = np.full_like(y_fit, np.nan, dtype=float)
    x_used[m] = apply_normalizer(x_fit[m], x_norm)
    y_used[m] = apply_normalizer(y_fit[m], y_norm)

    out[x_fit_name] = x_fit
    out[y_fit_name] = y_fit
    out["x_used"] = x_used
    out["y_used"] = y_used

    meta = {
        "space": space,
        "xcol_raw": xcol_raw,
        "ycol_raw": ycol_raw,
        "x_fit_name": x_fit_name,
        "y_fit_name": y_fit_name,
        "x_norm_method": x_norm["method"],
        "x_norm_center": x_norm["center"],
        "x_norm_scale": x_norm["scale"],
        "x_norm_vmin": x_norm.get("vmin", np.nan),
        "x_norm_vmax": x_norm.get("vmax", np.nan),
        "y_norm_method": y_norm["method"],
        "y_norm_center": y_norm["center"],
        "y_norm_scale": y_norm["scale"],
        "y_norm_vmin": y_norm.get("vmin", np.nan),
        "y_norm_vmax": y_norm.get("vmax", np.nan),
    }
    return out, meta

# =============================================================================
# MAIN
# =============================================================================

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("top", help="Top folder containing previous sip_pipeline outputs")
    ap.add_argument("--freq", type=float, default=0.01, help="Frequency tag")

    ap.add_argument("--include-log-space", action="store_true",
                    help="Also compare log(sig_real) vs log(S) in addition to real space.")

    ap.add_argument("--log-base", type=str, default="10", choices=["10", "e"],
                    help="Base used for transformed power/log-space work: 10 or natural log e.")

    ap.add_argument("--normalize", type=str, default="none",
                    choices=["none", "zscore", "minmax"],
                    help="Normalize x and y within each relevant group (sample/site/all).")

    ap.add_argument("--logistic-max-seconds", type=float, default=20.0,
                    help="Time cap per logistic fit")

    ap.add_argument("--min-points", type=int, default=4,
                    help="Minimum points required")

    ap.add_argument("--plot-equations", action="store_true",
                    help="Write equations and sigma/lnL/AIC on plots")

    ap.add_argument("--eq-position", type=str, default="right", choices=["right", "inside"],
                    help="Equation placement on plots")

    ap.add_argument("--no-plots", action="store_true",
                    help="Skip plot generation")

    ap.add_argument("--site-map", type=str, default="",
                    help="Optional CSV mapping run_name,site to override inferred site names")

    ap.add_argument("--outdir", type=str, default="COMPARE_COND_UNORM",
                    help="Output folder under top")

    args = ap.parse_args()

    top = Path(args.top).expanduser().resolve()
    if not top.exists():
        print(f"Not found: {top}")
        raise SystemExit(2)

    outdir = top / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    tag = tag_from_freq(args.freq)

    df_all = load_unorm_joined(top, tag)
    if df_all.empty:
        print(f"No joined__UNORM_{tag}Hz.csv files found.")
        raise SystemExit(1)

    if args.site_map:
        smap = pd.read_csv(args.site_map)
        if "run_name" in smap.columns and "site" in smap.columns:
            dmap = dict(zip(smap["run_name"].astype(str), smap["site"].astype(str)))
            df_all["__site"] = df_all["__run_name"].astype(str).map(
                lambda s: dmap.get(s, infer_site_name(s))
            )

    cases = list(DEFAULT_CASES)
    if args.include_log_space:
        cases.append(("LOG", "S", "sig_real"))

    stacked_csv = outdir / f"stacked_compare_input__UNORM_{tag}Hz.csv"
    df_all.to_csv(stacked_csv, index=False)

    compare_rows = []
    used_rows_all = []
    plot_jobs = []

    for space, xcol_raw, ycol_raw in cases:
        for level in ["SAMPLE", "SITE", "ALL"]:
            if level == "SAMPLE":
                groups = sorted(df_all["__run_name"].astype(str).unique().tolist())
                group_key = "__run_name"
            elif level == "SITE":
                groups = sorted(df_all["__site"].astype(str).unique().tolist())
                group_key = "__site"
            else:
                groups = ["ALL_DATA"]
                group_key = None

            for g in groups:
                if level == "ALL":
                    dfg = df_all.copy()
                else:
                    dfg = df_all.loc[df_all[group_key].astype(str) == g].copy()

                if dfg.empty:
                    continue

                prepared, meta = prepare_case_dataframe(
                    dfg=dfg,
                    space=space,
                    xcol_raw=xcol_raw,
                    ycol_raw=ycol_raw,
                    log_base=args.log_base,
                    normalize=args.normalize
                )
                if prepared is None:
                    continue

                x = safe_float_series(prepared["x_used"]).to_numpy(float)
                y = safe_float_series(prepared["y_used"]).to_numpy(float)
                m = np.isfinite(x) & np.isfinite(y)
                n_valid = int(np.sum(m))
                if n_valid < args.min_points:
                    continue

                fitcmp = compare_two_methods(
                    x=x[m],
                    y=y[m],
                    space=space,
                    log_base=args.log_base,
                    logistic_max_seconds=args.logistic_max_seconds,
                )
                if fitcmp is None:
                    continue

                folder_sources = "; ".join(sorted(prepared["__run_root"].astype(str).unique().tolist()))
                run_names = "; ".join(sorted(prepared["__run_name"].astype(str).unique().tolist()))

                p = fitcmp["power"]
                l = fitcmp["logistic4"]

                compare_rows.append({
                    "level": level,
                    "group": g,
                    "mode": "UNORM",
                    "space": space,
                    "xcol_raw": xcol_raw,
                    "ycol_raw": ycol_raw,
                    "n_valid": n_valid,
                    "winner": fitcmp["winner"],

                    "normalize": args.normalize,
                    "log_base": args.log_base,

                    "x_norm_method": meta["x_norm_method"],
                    "x_norm_center": meta["x_norm_center"],
                    "x_norm_scale": meta["x_norm_scale"],
                    "x_norm_vmin": meta["x_norm_vmin"],
                    "x_norm_vmax": meta["x_norm_vmax"],

                    "y_norm_method": meta["y_norm_method"],
                    "y_norm_center": meta["y_norm_center"],
                    "y_norm_scale": meta["y_norm_scale"],
                    "y_norm_vmin": meta["y_norm_vmin"],
                    "y_norm_vmax": meta["y_norm_vmax"],

                    "power_equation": p.get("equation", ""),
                    "power_a": p.get("a", np.nan),
                    "power_b": p.get("b", np.nan),
                    "power_r2": p.get("r2", np.nan),
                    "power_sigma_hat": p.get("sigma_hat", np.nan),
                    "power_lnL": p.get("lnL", np.nan),
                    "power_aic": p.get("aic", np.nan),
                    "power_sse": p.get("sse", np.nan),

                    "logistic_equation": l.get("equation", ""),
                    "logistic_amp_L": l.get("amp_L", np.nan),
                    "logistic_k": l.get("k", np.nan),
                    "logistic_x0": l.get("x0", np.nan),
                    "logistic_c": l.get("c", np.nan),
                    "logistic_r2": l.get("r2", np.nan),
                    "logistic_sigma_hat": l.get("sigma_hat", np.nan),
                    "logistic_lnL": l.get("lnL", np.nan),
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

                dfg_used = prepared.loc[m].copy()
                dfg_used["__compare_level"] = level
                dfg_used["__compare_group"] = g
                dfg_used["__space"] = space
                dfg_used["__winner"] = fitcmp["winner"]
                dfg_used["__normalize"] = args.normalize
                dfg_used["__log_base"] = args.log_base
                used_rows_all.append(dfg_used)

                if not args.no_plots:
                    png = outdir / (
                        f"{make_safe_name(level)}__{make_safe_name(g)}__"
                        f"{ycol_raw}_vs_{xcol_raw}__UNORM_{space}_{tag}Hz.png"
                    )
                    title = (
                        f"{level} {g}: {human_label(ycol_raw)} vs {human_label(xcol_raw)} "
                        f"[UNORM, {space}, {args.freq:g} Hz, norm={args.normalize}]"
                    )
                    plot_jobs.append((dfg_used.copy(), "x_used", "y_used", fitcmp, png, title))

    compare_df = pd.DataFrame(compare_rows)
    used_df = pd.concat(used_rows_all, axis=0, ignore_index=True) if used_rows_all else pd.DataFrame()

    for dfg_used, xcol_fit, ycol_fit, fitcmp, png, title in plot_jobs:
        try:
            save_compare_plot(
                dfg_used, xcol_fit, ycol_fit, fitcmp, png, title,
                plot_equations=args.plot_equations,
                eq_position=args.eq_position,
            )
        except Exception as e:
            print(f"[WARN] plot failed {png.name}: {e}", flush=True)

    compare_csv = outdir / f"compare_cond_unorm__{tag}Hz.csv"
    compare_df.to_csv(compare_csv, index=False)

    used_csv = outdir / f"compare_cond_unorm_used_data__{tag}Hz.csv"
    used_df.to_csv(used_csv, index=False)

    settings_df = pd.DataFrame({
        "key": [
            "top",
            "mode",
            "freq",
            "tag",
            "normalize",
            "log_base",
            "include_log_space",
            "logistic_max_seconds",
            "min_points",
            "plot_equations",
            "eq_position",
            "no_plots",
        ],
        "value": [
            str(top),
            "UNORM",
            args.freq,
            tag,
            args.normalize,
            args.log_base,
            args.include_log_space,
            args.logistic_max_seconds,
            args.min_points,
            args.plot_equations,
            args.eq_position,
            args.no_plots,
        ]
    })

    xlsx = outdir / f"compare_cond_unorm__{tag}Hz.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        if not used_df.empty:
            used_df.to_excel(w, sheet_name="used_data_ALL", index=False)
        if not compare_df.empty:
            compare_df.to_excel(w, sheet_name="compare_ALL", index=False)
        settings_df.to_excel(w, sheet_name="settings", index=False)

        for level in ["SAMPLE", "SITE", "ALL"]:
            dlev = compare_df.loc[compare_df["level"].astype(str) == level].copy() if not compare_df.empty else pd.DataFrame()
            if not dlev.empty:
                dlev.to_excel(w, sheet_name=sheet_safe_name(f"cmp_{level}"), index=False)

        for space in ["REAL", "LOG"]:
            dsp = compare_df.loc[compare_df["space"].astype(str) == space].copy() if not compare_df.empty else pd.DataFrame()
            if not dsp.empty:
                dsp.to_excel(w, sheet_name=sheet_safe_name(f"cmp_{space}"), index=False)

    print(f"Used data CSV: {used_csv}", flush=True)
    print(f"Compare CSV:   {compare_csv}", flush=True)
    print(f"Compare XLSX:  {xlsx}", flush=True)
    print(f"Plots folder:  {outdir}", flush=True)

if __name__ == "__main__":
    main()
