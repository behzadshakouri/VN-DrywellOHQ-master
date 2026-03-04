#!/usr/bin/env python3
"""
auto_sip_pipeline_v27_noscipy.py  (FULL — keep original plot style + fix logistic>=4 + no gnuplot warnings)

This version keeps the SAME plot STYLE as the first script you sent:
  - pngcairo sizes (1200x800 for single, 1200x1800 for diagnostics)
  - Arial,24
  - same labels/titles/keys positioning style
  - same gnuplot "column(\"name\")" pattern

Fixes included:
  1) Logistic4 can run with >=4 points (was >=5/6 before).
  2) No more gnuplot warnings / "x range is invalid" for logistic-only and diagnostics
     when logistic is not accepted:
       - We NEVER plot a "ps 0" dummy (gnuplot hates it).
       - We set fixed xrange/yrange and draw a simple y=0 line via inline data.
       - We only reference logistic columns in CSV when logistic is actually accepted.

Keeps:
  - both modes: NORM and UNORM (via --both-modes)
  - log-log + real-real fits (linear, poly, logistic4)
  - per-run all_values Excel + global all_values + summary CSV/XLSX
  - robust logbook matching by best overlap contiguous block
  - folder name as source label and S-1 -> S_{-1}

"""

import math
import argparse
import subprocess
import time
import re
from pathlib import Path

import numpy as np
import pandas as pd

# =============================================================================
# CONSTANTS / CONFIG
# =============================================================================

ANALYSIS_DIR = "Analysis"
RAW_DIR      = "SIP_Raw_Data"
RESULTS_DIR  = "Results"

ANALYSIS_XLSX_HINT = "output"
LOGBOOK_XLSX_HINT  = "Log Book"

LOGBOOK_ID_COLS_CAND     = ["Measurements","Measurement","Meas","MeasID","ID"]
LOGBOOK_SAT_COLS_CAND    = ["Degree Saturation M2","Degree Saturation","Saturation","Sw","S","DegreeSaturation"]
LOGBOOK_SAMPLE_COLS_CAND = ["Sample","Sample ID","SampleID","Soil Sample","Specimen"]

LOGBOOK_MASSWATER_CAND = ["Mass Water","Water Mass","Mass of Water","Mw","Water"]
LOGBOOK_POROSITY_CAND  = ["Soil Porosity","Porosity","phi","Φ"]
LOGBOOK_COLVOL_CAND    = ["Column Volume (m3)","Column Volume","Volume (m3)","Vcol"]

SIP_FREQ_CANDS = ["Frequency","Frequency (Hz)","Freq","Hz"]
SIP_RHO_CANDS  = ["Resistivity","Resistivity (Ω·m)","Resistivity (Ohm m)","rho"]
SIP_IMAG_CANDS = ["Imaginary Conductivity","Imaginary Conductivity (µS/cm)","Imag","sigma_imag"]

# Logistic bounds for stability
K_MIN = 0.05
K_MAX = 50.0

LOGISTIC_R2_EPS = 0.005
LOGISTIC_REQUIRE_BETTER_AIC = True

# =============================================================================
# BASIC HELPERS
# =============================================================================

def _norm(s: str) -> str:
    """Normalize column name comparisons: lowercase + strip."""
    return str(s).strip().lower()

def pick_col(df: pd.DataFrame, cands):
    """
    Return first matching dataframe column from candidates (case-insensitive).
    If none found, "".
    """
    low = {_norm(c): c for c in df.columns}
    for k in cands:
        kk = _norm(k)
        if kk in low:
            return low[kk]
    return ""

def find_first_matching_file(folder: Path, hint: str, ext=".xlsx"):
    """Find first file in folder tree with name containing hint (case-insensitive) and endswith ext."""
    if not folder.exists():
        return None
    for p in folder.rglob(f"*{ext}"):
        if hint.lower() in p.name.lower():
            return p
    return None

def to_float(x):
    """Robust float conversion (supports '40%' -> 40). Returns NaN on failure."""
    try:
        if isinstance(x, str):
            xs = x.strip()
            if xs.endswith("%"):
                return float(xs[:-1].strip())
            return float(xs)
        return float(x)
    except Exception:
        return float("nan")

def safe_log10(x):
    """Safe log10; returns NaN if x not finite or <=0."""
    x = float(x)
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return math.log10(x)

def sse(y, yhat):
    """Sum of squared errors."""
    r = y - yhat
    return float(np.sum(r*r))

def _is_blank(x) -> bool:
    """True if value is blank-ish (None, '', 'nan', 'none')."""
    if x is None:
        return True
    s = str(x).strip()
    return (s == "") or (s.lower() in ["nan", "none"])

def norm_measid(s: str) -> str:
    """
    Normalize measurement IDs so formatting typos don't break matching.
      "M 1" -> "m1"
      "M-1" -> "m1"
      "M_1" -> "m1"
    """
    s = str(s).strip().lower()
    s = re.sub(r'[^a-z0-9]+', '', s)
    return s

def pretty_source_name(run_folder_name: str) -> str:
    """
    Readable plot title from folder name.
    Also convert S-1 -> S_{-1} for gnuplot enhanced text.
    """
    s = str(run_folder_name).replace("_", " ").strip()
    s = re.sub(r'\bS-(\d+)\b', r'S_{-\1}', s)
    return s

def tag_from_freq(freq: float) -> str:
    """0.01 -> 0p01 (safe in filenames)."""
    return str(freq).replace(".", "p")

# =============================================================================
# F-DISTRIBUTION + P-VALUE (NO SCIPY)
# =============================================================================

def betacf(a, b, x, maxit=200, eps=3e-14):
    """Continued fraction approximation for incomplete beta."""
    am=1.0; bm=1.0; az=1.0
    qab=a+b; qap=a+1.0; qam=a-1.0
    bz=1.0 - qab*x/qap
    for m in range(1, maxit+1):
        em=float(m); tem=em+em
        d = em*(b-em)*x/((qam+tem)*(a+tem))
        ap = az + d*am
        bp = bz + d*bm
        d = -(a+em)*(qab+em)*x/((a+tem)*(qap+tem))
        app = ap + d*az
        bpp = bp + d*bz
        aold=az
        am=ap/bpp
        bm=bp/bpp
        az=app/bpp
        bz=1.0
        if abs(az-aold) < eps*abs(az):
            return az
    return az

def betai(a, b, x):
    """Regularized incomplete beta I_x(a,b)."""
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    ln_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a+b)
    bt = math.exp(math.log(x)*a + math.log(1.0-x)*b - ln_beta)
    if x < (a+1.0)/(a+b+2.0):
        return bt*betacf(a,b,x)/a
    else:
        return 1.0 - bt*betacf(b,a,1.0-x)/b

def f_dist_sf(F, d1, d2):
    """Survival function (1-CDF) for F distribution."""
    if not np.isfinite(F) or F < 0:
        return np.nan
    x = (d1*F)/(d1*F + d2)
    cdf = betai(d1/2.0, d2/2.0, x)
    return max(0.0, 1.0 - cdf)

def f_test_pvalue(ss_res, ss_tot, df_model, n):
    """
    Model p-value via F-test.
    df_model = predictors count excluding intercept.
      linear y=a*x+b -> df_model=1
      log4 approx -> df_model=3
    """
    if not (np.isfinite(ss_res) and np.isfinite(ss_tot)):
        return np.nan
    if ss_tot <= 0:
        return np.nan
    if ss_res >= ss_tot:
        return 1.0
    df2 = n - df_model - 1
    if df2 <= 0:
        return np.nan
    num = (ss_tot - ss_res)/df_model
    den = ss_res/df2
    if den <= 0:
        return np.nan
    F = num/den
    return f_dist_sf(F, df_model, df2)

def aic_gaussian(ss_res, n, k_params):
    """Gaussian AIC: AIC = n*ln(ss_res/n) + 2*k_params."""
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)

# =============================================================================
# FIT ROUTINES: LINEAR / THROUGH ORIGIN / POLY
# =============================================================================

def linear_fit(x, y):
    """
    OLS y = slope*x + intercept.
    Returns slope, intercept, r2, pvalue, ss_res, aic
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    yhat = slope*x + intercept
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    pval = f_test_pvalue(ss_res, ss_tot, df_model=1, n=len(x))
    aic = aic_gaussian(ss_res, len(x), k_params=2)
    return float(slope), float(intercept), float(r2), float(pval), float(ss_res), float(aic)

def linear_fit_through_origin(x, y):
    """Least squares forced through origin: y=a*x."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 2:
        return np.nan, np.nan, np.nan
    den = float(np.dot(x, x))
    if den <= 0:
        return np.nan, np.nan, np.nan
    a = float(np.dot(x, y) / den)
    yhat = a*x
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    return a, float(r2), float(ss_res)

def poly_fit_predict(x, y, deg: int):
    """
    Polynomial regression y ~ poly(x,deg).
    Returns (coeffs_high_to_low, yhat, r2, ss_res).
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < (deg + 1):
        return None, None, np.nan, np.nan
    try:
        coeffs = np.polyfit(x, y, deg)
        yhat = np.polyval(coeffs, x)
        ss_res = sse(y, yhat)
        ss_tot = sse(y, np.full_like(y, np.mean(y)))
        r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
        return coeffs, yhat, float(r2), float(ss_res)
    except Exception:
        return None, None, np.nan, np.nan

# =============================================================================
# LOGISTIC FIT (NO SCIPY) - 4 PARAM
# =============================================================================

def logistic4(x, L, k, x0, c):
    """f(x) = c + L/(1+exp(-k*(x-x0)))."""
    z = k*(x - x0)
    z = np.clip(z, -60.0, 60.0)
    return c + L/(1.0 + np.exp(-z))

def _logistic_init_guess(x, y):
    """Heuristic init for log4."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    idx = np.argsort(x)
    x = x[idx]; y = y[idx]

    y_min = float(np.min(y))
    y_max = float(np.max(y))
    L = max(1e-6, y_max - y_min)
    c = y_min

    y_mid = y_min + 0.5*L
    j = int(np.argmin(np.abs(y - y_mid)))
    x0 = float(x[j])

    if len(x) >= 3:
        j0 = max(1, min(len(x)-2, j))
        dy = float(y[j0+1] - y[j0-1])
        dx = float(x[j0+1] - x[j0-1])
        slope = dy/dx if dx != 0 else 0.0
        k = abs(4.0*slope/max(L,1e-6))
        k = float(np.clip(k, K_MIN, K_MAX))
    else:
        k = 1.0

    return L, k, x0, c

def logistic_fit_noscipy(x, y, fast=False, max_seconds=0.0):
    """
    Coarse coordinate search around init guess (SciPy-free).
    UPDATED: allow >=4 points (your request).
    Returns (fit_dict, status_str).
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    n = len(x)

    if n < 4:
        return None, "not_enough_points"
    if not fast:
        return None, "disabled"

    t0 = time.time()
    Lb, kb, x0b, cb = _logistic_init_guess(x, y)

    x_min, x_max = float(np.min(x)), float(np.max(x))
    spanx = max(1e-6, x_max - x_min)

    stepL  = 0.25*max(abs(Lb), 1e-6)
    stepk  = 0.25*max(abs(kb), 1e-6)
    stepx0 = 0.20*spanx
    stepc  = 0.25*max(abs(Lb), 1e-6)

    def clamp_params(L,k,x0,c):
        L = max(1e-9, L)
        k = float(np.clip(k, K_MIN, K_MAX))
        return L,k,x0,c

    def current_sse(L,k,x0,c):
        return sse(y, logistic4(x, L,k,x0,c))

    Lb,kb,x0b,cb = clamp_params(Lb,kb,x0b,cb)
    best = current_sse(Lb,kb,x0b,cb)

    for _ in range(10):
        if max_seconds and (time.time() - t0) > max_seconds:
            return None, "timeout"
        improved = False
        for dL in [0, -stepL, stepL]:
            for dk in [0, -stepk, stepk]:
                for dx0 in [0, -stepx0, stepx0]:
                    for dc in [0, -stepc, stepc]:
                        if dL==dk==dx0==dc==0:
                            continue
                        if max_seconds and (time.time() - t0) > max_seconds:
                            return None, "timeout"
                        L2,k2,x02,c2 = clamp_params(Lb + dL, kb + dk, x0b + dx0, cb + dc)
                        s2 = current_sse(L2,k2,x02,c2)
                        if s2 + 1e-12 < best:
                            best = s2
                            Lb,kb,x0b,cb = L2,k2,x02,c2
                            improved = True
        stepL *= 0.6; stepk *= 0.6; stepx0 *= 0.6; stepc *= 0.6
        if not improved and max(stepL, stepk, stepx0, stepc) < 1e-4:
            break

    yhat = logistic4(x, Lb, kb, x0b, cb)
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    p_model = f_test_pvalue(ss_res, ss_tot, df_model=3, n=n)
    aic = aic_gaussian(ss_res, n, k_params=4)

    return {
        "L": float(Lb), "k": float(kb), "x0": float(x0b), "c": float(cb),
        "r2": float(r2),
        "p_model": float(p_model) if np.isfinite(p_model) else p_model,
        "aic": float(aic) if np.isfinite(aic) else aic,
        "ss_res": float(ss_res), "ss_tot": float(ss_tot)
    }, "ok"

def accept_logistic(lin_r2, lin_aic, log_fit):
    """Accept logic: R2 improves or AIC improves."""
    if log_fit is None:
        return False, "no_fit"
    log_r2 = log_fit.get("r2", np.nan)
    log_aic = log_fit.get("aic", np.nan)
    if not np.isfinite(log_r2) or not np.isfinite(log_aic):
        return False, "bad_metrics"
    if np.isfinite(lin_r2) and (log_r2 >= lin_r2 + LOGISTIC_R2_EPS):
        return True, "R2_improved"
    if LOGISTIC_REQUIRE_BETTER_AIC and np.isfinite(lin_aic) and (log_aic < lin_aic):
        return True, "AIC_improved"
    return False, "rejected"

# =============================================================================
# SIP OUTPUT READING
# =============================================================================

def read_sip_output(sip_xlsx: Path, target_freq: float, max_df: float, sheet_filter: str|None):
    """
    Read SIP output Excel:
      each sheet -> measurement
      pick row closest to target_freq
      return dataframe with measurement, measurement_norm, freq_hz, rho, sigma_imag
    """
    xls = pd.ExcelFile(sip_xlsx)
    out = []
    for sh in xls.sheet_names:
        if sheet_filter and (sheet_filter.lower() not in str(sh).lower()):
            continue
        df = pd.read_excel(sip_xlsx, sheet_name=sh)
        if df is None or df.empty:
            continue

        freq_col = pick_col(df, SIP_FREQ_CANDS) or df.columns[0]
        rho_col  = pick_col(df, SIP_RHO_CANDS)
        imag_col = pick_col(df, SIP_IMAG_CANDS)

        # fallback indices if names not found
        if (not rho_col) or (not imag_col):
            if len(df.columns) >= 5:
                imag_col = df.columns[1]
                rho_col  = df.columns[4]
            else:
                continue

        freqs = df[freq_col].apply(to_float).to_numpy(float)
        if not np.isfinite(freqs).any():
            continue

        idx = int(np.nanargmin(np.abs(freqs - target_freq)))
        f_found = float(freqs[idx])
        if max_df > 0 and abs(f_found - target_freq) > max_df:
            continue

        rho = to_float(df.iloc[idx][rho_col])
        sim = to_float(df.iloc[idx][imag_col])

        out.append({
            "measurement": str(sh).strip(),
            "measurement_norm": norm_measid(sh),
            "freq_hz": f_found,
            "rho": float(rho),
            "sigma_imag": float(sim)
        })

    return pd.DataFrame(out)

# =============================================================================
# LOGBOOK READING: BEST BLOCK BY OVERLAP
# =============================================================================

def _contiguous_blocks(mask: np.ndarray):
    """Return contiguous (i0,i1) blocks where mask True."""
    blocks = []
    in_block = False
    i0 = 0
    for i, v in enumerate(mask):
        if v and (not in_block):
            in_block = True
            i0 = i
        elif (not v) and in_block:
            blocks.append((i0, i-1))
            in_block = False
    if in_block:
        blocks.append((i0, len(mask)-1))
    return blocks

def read_logbook_bestblock_by_overlap(log_xlsx: Path, sip_measurements_norm: list[str]):
    """
    Choose best logbook sheet/block maximizing overlap with SIP measurement IDs.
    Returns sat_map, used_sheet, overlap_cnt, block_rows.
    """
    xls = pd.ExcelFile(log_xlsx)
    sip_set = set([m for m in sip_measurements_norm if not _is_blank(m)])

    best = None
    for sh in xls.sheet_names:
        try:
            df = pd.read_excel(log_xlsx, sheet_name=sh)
        except Exception:
            continue
        if df is None or df.empty:
            continue

        id_col  = pick_col(df, LOGBOOK_ID_COLS_CAND) or ""
        sat_col = pick_col(df, LOGBOOK_SAT_COLS_CAND)
        mw_col  = pick_col(df, LOGBOOK_MASSWATER_CAND)
        por_col = pick_col(df, LOGBOOK_POROSITY_CAND)
        vol_col = pick_col(df, LOGBOOK_COLVOL_CAND)

        if not id_col:
            continue

        has_meas = ~df[id_col].astype(str).apply(_is_blank)

        if sat_col:
            has_sat = df[sat_col].apply(lambda x: np.isfinite(to_float(x)))
        else:
            if not (mw_col and por_col and vol_col):
                continue
            has_sat = df.apply(
                lambda r: np.isfinite(to_float(r.get(mw_col, np.nan)))
                          and np.isfinite(to_float(r.get(por_col, np.nan)))
                          and np.isfinite(to_float(r.get(vol_col, np.nan))),
                axis=1
            )

        usable = (has_meas & has_sat).to_numpy(bool)
        blocks = _contiguous_blocks(usable)
        if not blocks:
            continue

        for (i0, i1) in blocks:
            blk = df.iloc[i0:i1+1].copy()
            blk_ids = blk[id_col].astype(str).apply(norm_measid).tolist()
            overlap = len(set(blk_ids) & sip_set)
            score = overlap * 1000 + (i1 - i0 + 1)
            if (best is None) or (score > best[0]):
                best = (score, sh, blk, id_col, sat_col, mw_col, por_col, vol_col, overlap)

    if best is None:
        raise RuntimeError(f"Logbook has no usable saturation blocks. Sheets: {xls.sheet_names}")

    _, sh, blk, id_col, sat_col, mw_col, por_col, vol_col, overlap = best

    sat_map = {}
    for _, r in blk.iterrows():
        mid = str(r.get(id_col, "")).strip()
        if _is_blank(mid):
            continue

        S = np.nan
        if sat_col:
            S = to_float(r.get(sat_col, np.nan))
            if np.isfinite(S) and S > 1.0 and S <= 100.0:
                S = S/100.0
        else:
            Mw  = to_float(r.get(mw_col, np.nan))   # g
            phi = to_float(r.get(por_col, np.nan))
            V   = to_float(r.get(vol_col, np.nan))  # m3
            if np.isfinite(Mw) and np.isfinite(phi) and np.isfinite(V) and phi > 0 and V > 0:
                Mw_kg = Mw/1000.0
                rho_w = 1000.0
                S = Mw_kg / (rho_w * (phi * V))

        if np.isfinite(S) and S > 0:
            sat_map[norm_measid(mid)] = float(S)

    if not sat_map:
        raise RuntimeError(f"Selected logbook sheet '{sh}' but mapped 0 saturation rows.")

    return sat_map, sh, int(overlap), int(len(blk))

# =============================================================================
# RESULTS CLEANING
# =============================================================================

def clean_results_dir(results_dir: Path):
    """Delete files in Results/ (keep folder)."""
    if not results_dir.exists():
        return
    for p in results_dir.iterdir():
        try:
            if p.is_file() or p.is_symlink():
                p.unlink()
        except Exception:
            pass

# =============================================================================
# GNUPLOT EXEC
# =============================================================================

def _run_gnuplot(script_name: str, cwd: Path, timeout_s: int, fatal: bool):
    """Run gnuplot script in cwd."""
    try:
        subprocess.run(["gnuplot", script_name], cwd=str(cwd), check=True,
                       timeout=(timeout_s if timeout_s > 0 else None))
        return True, ""
    except Exception as e:
        msg = str(e)
        if fatal:
            raise
        print(f"[WARN] gnuplot failed for {cwd}/{script_name}: {msg}", flush=True)
        return False, msg

# =============================================================================
# GNUPLOT WRITERS (KEEP STYLE)
# =============================================================================

def write_plot_logistic_only(results_dir: Path, csv_name: str, out_png: str,
                             target_freq: float, source_label: str,
                             logistic_accepted: bool):
    """
    SAME STYLE as your original; only fixes the dummy plot so gnuplot never errors.
    """
    label = "{:g}".format(target_freq)

    if logistic_accepted:
        gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,24'
set output '{out_png}'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
set title '{source_label}: logistic-only (accepted) @ {label} Hz'
plot '{csv_name}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logRho_logfit")) with lines dt 1 lw 4 title 'logistic (accepted)'
"""
    else:
        # IMPORTANT: draw something valid (no ps 0), and set explicit ranges.
        gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,24'
set output '{out_png}'
set grid
unset key
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
set title '{source_label}: logistic-only @ {label} Hz'
set label 1 'No logistic accepted' at graph 0.5,0.85 center front
set xrange [-2:0]
set yrange [*:*]
plot '-' using 1:2 with lines lw 3 notitle
-2 0
0  0
e
"""
    (results_dir / "plot_logistic_only.gp").write_text(gp)

def write_plot_logistic_diag(results_dir: Path, csv_name: str, out_png: str,
                             target_freq: float, source_label: str,
                             logistic_accepted: bool):
    """
    SAME STYLE as your original diagnostics plot; only fix:
      - no CSV logistic columns referenced unless logistic_accepted
      - dummy panel uses a real line + fixed ranges (no ps 0)
    """
    label = "{:g}".format(target_freq)

    extra_curve = ""
    if logistic_accepted:
        extra_curve = (
            ", \\\n     '{csv}' using (column(\"logS\")):(column(\"logRho_logfit\")) "
            "with lines dt 1 lw 4 title 'logistic (accepted)'"
        ).format(csv=csv_name)

    if logistic_accepted:
        extra3 = ""
        panel3 = (
            "plot '{csv}' using (column(\"logS\")):(column(\"resid_log_rho\")) "
            "with points pt 7 ps 1.4 title 'residuals'\n"
        ).format(csv=csv_name)
    else:
        extra3 = (
            "unset key\n"
            "set label 1 'No logistic accepted' at graph 0.5,0.5 center front\n"
            "set xrange [-2:0]\n"
            "set yrange [-1:1]\n"
        )
        panel3 = (
            "plot '-' using 1:2 with lines lw 3 notitle\n"
            "-2 0\n"
            "0  0\n"
            "e\n"
        )

    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1800 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 3,1

set title '{source_label}: Diagnostics: logRho vs logS at {label} Hz'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
plot '{csv_name}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logRho_linfit")) with lines dt 2 lw 3 title 'linear'{extra_curve}

set title 'Residuals: linear (data - linear)'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'residual'
plot '{csv_name}' using (column("logS")):(column("resid_lin_rho")) with points pt 7 ps 1.4 title 'residuals'

set title 'Residuals: logistic (data - logistic)'
set grid
set xlabel 'log10(Saturation)'
set ylabel 'residual'
{extra3}
{panel3}

unset multiplot
"""
    (results_dir / "plot_logistic_diagnostics.gp").write_text(gp)

def write_plot_xy_bundle(results_dir: Path, csv_name: str, out_png: str, out_gp: str,
                         source_label: str, freq: float,
                         x_expr: str, y_expr: str,
                         xlab: str, ylab: str,
                         title: str,
                         do_linear: bool, do_poly: bool, poly_deg: int,
                         do_log4: bool,
                         y_lin_col: str, y_poly_col: str, y_log4_col: str):
    """
    Same as your original generator: precomputed columns, column("name") gnuplot.
    Avoids plotting missing curves by checking if column has finite data.
    """
    label = "{:g}".format(freq)

    df = None
    try:
        df = pd.read_csv(results_dir / csv_name)
    except Exception:
        df = None

    def col_has_finite(colname: str) -> bool:
        if df is None: return True
        if colname not in df.columns: return False
        v = pd.to_numeric(df[colname], errors="coerce").to_numpy(float)
        return bool(np.isfinite(v).any())

    lines = []
    lines.append(f"'{csv_name}' using ({x_expr}):({y_expr}) with points pt 7 ps 1.7 title 'data'")
    if do_linear and y_lin_col and col_has_finite(y_lin_col):
        lines.append(f"'{csv_name}' using ({x_expr}):(column(\"{y_lin_col}\")) with lines dt 2 lw 3 title 'linear'")
    if do_poly and poly_deg >= 2 and y_poly_col and col_has_finite(y_poly_col):
        lines.append(f"'{csv_name}' using ({x_expr}):(column(\"{y_poly_col}\")) with lines dt 1 lw 4 title 'poly deg {poly_deg}'")
    if do_log4 and y_log4_col and col_has_finite(y_log4_col):
        lines.append(f"'{csv_name}' using ({x_expr}):(column(\"{y_log4_col}\")) with lines dt 3 lw 4 title 'logistic4'")

    plot_cmd = "plot " + ", \\\n     ".join(lines) + "\n"

    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,24'
set output '{out_png}'
set grid
set key top right
set tics out

set title '{title} ({label} Hz)'
set xlabel '{xlab}'
set ylabel '{ylab}'

{plot_cmd}
"""
    (results_dir / out_gp).write_text(gp)

# =============================================================================
# DISCOVERY
# =============================================================================

def discover_run_roots_fast(base: Path):
    """Run roots are folders containing Analysis/ and SIP_Raw_Data/."""
    roots = set()
    for a in base.rglob(ANALYSIS_DIR):
        if not a.is_dir():
            continue
        run = a.parent
        if (run / RAW_DIR).is_dir():
            roots.add(run)
    return sorted(roots)

# =============================================================================
# MAIN PROCESSING FOR ONE RUN (ONE MODE)
# =============================================================================

def process_one_mode(run_root: Path, target_freq: float, sheet_filter: str|None,
                     min_points: int, max_df: float,
                     do_log4: bool, logistic_fast: bool, logistic_max_seconds: float,
                     gnuplot_timeout: int,
                     clean_results: bool,
                     ici_mode: str,
                     poly_deg: int,
                     poly_raw: bool,
                     mode: str,
                     all_rows_accum: list,
                     per_run_excel_rows: list):
    """
    Process one run folder for one mode:
      mode = "NORM" or "UNORM"
    Produces joined__<mode>_<tag>Hz.csv and the plots for this mode.
    Also writes legacy joined_<tag>Hz.csv for NORM compatibility.
    """
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return {"run": str(run_root), "mode": mode, "status": "SKIP", "reason": "missing excel", "points": 0}

    source_label = pretty_source_name(run_root.name)
    results_dir.mkdir(parents=True, exist_ok=True)
    if clean_results and mode == "NORM":
        clean_results_dir(results_dir)

    tag = tag_from_freq(target_freq)

    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df, sheet_filter=sheet_filter)
    if sip_df.empty:
        return {"run": str(run_root), "mode": mode, "status": "FAIL", "reason": "no sip sheets", "points": 0}

    sat_map, used_sheet, overlap_cnt, block_rows = read_logbook_bestblock_by_overlap(
        log_xlsx, sip_df["measurement_norm"].tolist()
    )
    sip_df["S"] = sip_df["measurement_norm"].map(sat_map).astype(float)

    sip_df["logS"]    = sip_df["S"].apply(safe_log10)
    sip_df["logRho"]  = sip_df["rho"].apply(safe_log10)
    sip_df["logSig"]  = sip_df["sigma_imag"].apply(safe_log10)

    dbg = sip_df.copy()
    dbg.to_csv(results_dir / f"debug_before_filter_{tag}Hz__{mode}.csv", index=False)

    mask = np.isfinite(sip_df["logS"]) & np.isfinite(sip_df["logRho"]) & np.isfinite(sip_df["logSig"])
    dropped = sip_df.loc[~mask, ["measurement", "S", "rho", "sigma_imag"]].copy()
    if len(dropped) > 0:
        dropped.to_csv(results_dir / f"dropped_rows_{tag}Hz__{mode}.csv", index=False)
    sip_df = sip_df.loc[mask].copy()

    if len(sip_df) < min_points:
        return {"run": str(run_root), "mode": mode, "status": "FAIL", "reason": f"not enough points ({len(sip_df)})", "points": int(len(sip_df))}

    iref = int(np.nanargmax(sip_df["S"].to_numpy(float)))
    rho_ref  = float(sip_df.iloc[iref]["rho"])
    sig_ref  = float(sip_df.iloc[iref]["sigma_imag"])
    logRho_ref = safe_log10(rho_ref)
    logSig_ref = safe_log10(sig_ref)

    sip_df["rho_ref"] = rho_ref
    sip_df["sig_ref"] = sig_ref
    sip_df["logRho_ref"] = logRho_ref
    sip_df["logSig_ref"] = logSig_ref

    if mode == "NORM":
        sip_df["rho_idx"]   = sip_df["rho"] / rho_ref
        sip_df["sig_idx"]   = sip_df["sigma_imag"] / sig_ref
        sip_df["logRhoIdx"] = sip_df["logRho"] - logRho_ref

        if ici_mode == "sig_over_sigref":
            sip_df["logSigIdx"] = sip_df["logSig"] - logSig_ref
        else:
            sip_df["logSigIdx"] = logSig_ref - sip_df["logSig"]

        sip_df["logRI_vn"]  = sip_df["logRhoIdx"]
        sip_df["logICI_vn"] = logSig_ref - sip_df["logSig"]
    else:
        sip_df["rho_idx"]   = sip_df["rho"]
        sip_df["sig_idx"]   = sip_df["sigma_imag"]
        sip_df["logRhoIdx"] = sip_df["logRho"]
        sip_df["logSigIdx"] = sip_df["logSig"]
        sip_df["logRI_vn"]  = sip_df["logRho"]
        sip_df["logICI_vn"] = sip_df["logSig"]

    # =========================
    # LOG-LOG FITS (x=logS)
    # =========================
    a_lr, b_lr, r2_lr, p_lr, ss_lr, aic_lr = linear_fit(sip_df["logS"], sip_df["logRho"])
    a_ls, b_ls, r2_ls, p_ls, ss_ls, aic_ls = linear_fit(sip_df["logS"], sip_df["logSig"])

    sip_df["logRho_lin"] = a_lr*sip_df["logS"] + b_lr
    sip_df["logSig_lin"] = a_ls*sip_df["logS"] + b_ls

    sip_df["logRho_poly"] = np.nan
    sip_df["logSig_poly"] = np.nan
    r2_logrho_poly = np.nan
    r2_logsig_poly = np.nan
    if poly_deg and poly_deg >= 2:
        cR, _, r2pR, _ = poly_fit_predict(sip_df["logS"], sip_df["logRho"], poly_deg)
        cS, _, r2pS, _ = poly_fit_predict(sip_df["logS"], sip_df["logSig"], poly_deg)
        r2_logrho_poly = r2pR
        r2_logsig_poly = r2pS
        if cR is not None:
            sip_df["logRho_poly"] = np.polyval(cR, sip_df["logS"].to_numpy(float))
        if cS is not None:
            sip_df["logSig_poly"] = np.polyval(cS, sip_df["logS"].to_numpy(float))

    sip_df["logRho_log4"] = np.nan
    sip_df["logSig_log4"] = np.nan
    log4R = None
    log4S = None

    logistic_attempted = False
    logistic_accepted = False
    logistic_reason = ""
    logistic_status = ""

    if do_log4 and len(sip_df) >= 4:
        logistic_attempted = True
        log4R, stR = logistic_fit_noscipy(sip_df["logS"], sip_df["logRho"], fast=logistic_fast, max_seconds=logistic_max_seconds)
        logistic_status = stR
        if log4R is not None:
            sip_df["logRho_log4"] = logistic4(sip_df["logS"].to_numpy(float), log4R["L"], log4R["k"], log4R["x0"], log4R["c"])

        ok, why = accept_logistic(r2_lr, aic_lr, log4R)
        logistic_accepted = bool(ok)
        logistic_reason = why

        log4S, _ = logistic_fit_noscipy(sip_df["logS"], sip_df["logSig"], fast=logistic_fast, max_seconds=logistic_max_seconds)
        if log4S is not None:
            sip_df["logSig_log4"] = logistic4(sip_df["logS"].to_numpy(float), log4S["L"], log4S["k"], log4S["x0"], log4S["c"])
    else:
        logistic_attempted = bool(do_log4)
        logistic_status = "not_enough_points" if (do_log4 and len(sip_df) < 4) else "disabled"

    # residuals for diagnostics
    sip_df["logRho_linfit"] = sip_df["logRho_lin"]
    sip_df["resid_lin_rho"] = sip_df["logRho"] - sip_df["logRho_linfit"]

    if logistic_accepted:
        sip_df["logRho_logfit"] = sip_df["logRho_log4"]
        sip_df["resid_log_rho"] = sip_df["logRho"] - sip_df["logRho_logfit"]
    else:
        sip_df["logRho_logfit"] = np.nan
        sip_df["resid_log_rho"] = np.nan

    # =========================
    # REAL-REAL FITS (x=S)
    # =========================
    sip_df["S_real"] = sip_df["S"].astype(float)

    if mode == "NORM":
        y_rho_real = sip_df["rho_idx"].to_numpy(float)
        y_sig_real = sip_df["sig_idx"].to_numpy(float)
    else:
        y_rho_real = sip_df["rho"].to_numpy(float)
        y_sig_real = sip_df["sigma_imag"].to_numpy(float)

    xS = sip_df["S_real"].to_numpy(float)

    mr, cr, r2r, pr, ssr, aicr = linear_fit(xS, y_rho_real)
    ms, cs, r2s, ps, sss, aics = linear_fit(xS, y_sig_real)

    sip_df["rho_real"] = y_rho_real
    sip_df["sig_real"] = y_sig_real
    sip_df["rho_lin_real"] = mr*sip_df["S_real"] + cr
    sip_df["sig_lin_real"] = ms*sip_df["S_real"] + cs

    sip_df["rho_poly_real"] = np.nan
    sip_df["sig_poly_real"] = np.nan
    r2_rho_poly_real = np.nan
    r2_sig_poly_real = np.nan
    if poly_deg and poly_deg >= 2:
        cpr, _, r2pr, _ = poly_fit_predict(xS, y_rho_real, poly_deg)
        cps, _, r2ps, _ = poly_fit_predict(xS, y_sig_real, poly_deg)
        r2_rho_poly_real = r2pr
        r2_sig_poly_real = r2ps
        if cpr is not None:
            sip_df["rho_poly_real"] = np.polyval(cpr, xS)
        if cps is not None:
            sip_df["sig_poly_real"] = np.polyval(cps, xS)

    sip_df["rho_log4_real"] = np.nan
    sip_df["sig_log4_real"] = np.nan
    if do_log4 and len(sip_df) >= 4:
        log4R_real, _ = logistic_fit_noscipy(xS, y_rho_real, fast=logistic_fast, max_seconds=logistic_max_seconds)
        log4S_real, _ = logistic_fit_noscipy(xS, y_sig_real, fast=logistic_fast, max_seconds=logistic_max_seconds)
        if log4R_real is not None:
            sip_df["rho_log4_real"] = logistic4(xS, log4R_real["L"], log4R_real["k"], log4R_real["x0"], log4R_real["c"])
        if log4S_real is not None:
            sip_df["sig_log4_real"] = logistic4(xS, log4S_real["L"], log4S_real["k"], log4S_real["x0"], log4S_real["c"])

    # =========================
    # WRITE JOINED CSV
    # =========================
    joined_name = f"joined__{mode}_{tag}Hz.csv"
    joined_path = results_dir / joined_name
    sip_df.to_csv(joined_path, index=False)

    # Legacy compatibility: keep old joined_<tag>Hz.csv for NORM
    if mode == "NORM":
        (results_dir / f"joined_{tag}Hz.csv").write_text(joined_path.read_text())

    # =========================
    # PLOTS: REAL + LOGLOG for rho and conductivity
    # =========================
    methods = "LINEAR"
    if poly_deg and poly_deg >= 2:
        methods += "+POLY"
    if do_log4:
        methods += "+LOG4"

    # LOGLOG
    write_plot_xy_bundle(
        results_dir=results_dir,
        csv_name=joined_name,
        out_png=f"rho_vs_S__LOGLOG__{mode}__{methods}_{tag}Hz.png",
        out_gp =f"rho_vs_S__LOGLOG__{mode}__{methods}_{tag}Hz.gp",
        source_label=source_label,
        freq=target_freq,
        x_expr='column("logS")',
        y_expr='column("logRho")',
        xlab='log10(Saturation)',
        ylab='log10(Resistivity)',
        title=f"{source_label}: rho vs S (LOGLOG) [{mode}]",
        do_linear=True,
        do_poly=(poly_deg and poly_deg >= 2),
        poly_deg=poly_deg,
        do_log4=do_log4,
        y_lin_col="logRho_lin",
        y_poly_col="logRho_poly",
        y_log4_col="logRho_log4"
    )

    write_plot_xy_bundle(
        results_dir=results_dir,
        csv_name=joined_name,
        out_png=f"cond_vs_S__LOGLOG__{mode}__{methods}_{tag}Hz.png",
        out_gp =f"cond_vs_S__LOGLOG__{mode}__{methods}_{tag}Hz.gp",
        source_label=source_label,
        freq=target_freq,
        x_expr='column("logS")',
        y_expr='column("logSig")',
        xlab='log10(Saturation)',
        ylab='log10(Conductivity)',
        title=f"{source_label}: conductivity vs S (LOGLOG) [{mode}]",
        do_linear=True,
        do_poly=(poly_deg and poly_deg >= 2),
        poly_deg=poly_deg,
        do_log4=do_log4,
        y_lin_col="logSig_lin",
        y_poly_col="logSig_poly",
        y_log4_col="logSig_log4"
    )

    # REAL
    write_plot_xy_bundle(
        results_dir=results_dir,
        csv_name=joined_name,
        out_png=f"rho_vs_S__REAL__{mode}__{methods}_{tag}Hz.png",
        out_gp =f"rho_vs_S__REAL__{mode}__{methods}_{tag}Hz.gp",
        source_label=source_label,
        freq=target_freq,
        x_expr='column("S_real")',
        y_expr='column("rho_real")',
        xlab='Saturation',
        ylab=('rho/rho_ref' if mode == "NORM" else 'Resistivity'),
        title=f"{source_label}: rho vs S (REAL) [{mode}]",
        do_linear=True,
        do_poly=(poly_deg and poly_deg >= 2),
        poly_deg=poly_deg,
        do_log4=do_log4,
        y_lin_col="rho_lin_real",
        y_poly_col="rho_poly_real",
        y_log4_col="rho_log4_real"
    )

    write_plot_xy_bundle(
        results_dir=results_dir,
        csv_name=joined_name,
        out_png=f"cond_vs_S__REAL__{mode}__{methods}_{tag}Hz.png",
        out_gp =f"cond_vs_S__REAL__{mode}__{methods}_{tag}Hz.gp",
        source_label=source_label,
        freq=target_freq,
        x_expr='column("S_real")',
        y_expr='column("sig_real")',
        xlab='Saturation',
        ylab=('sig/sig_ref' if mode == "NORM" else 'Conductivity'),
        title=f"{source_label}: conductivity vs S (REAL) [{mode}]",
        do_linear=True,
        do_poly=(poly_deg and poly_deg >= 2),
        poly_deg=poly_deg,
        do_log4=do_log4,
        y_lin_col="sig_lin_real",
        y_poly_col="sig_poly_real",
        y_log4_col="sig_log4_real"
    )

    # Diagnostics + logistic-only (legacy names; NORM only)
    if mode == "NORM":
        legacy_csv = f"joined_{tag}Hz.csv"
        write_plot_logistic_only(
            results_dir, legacy_csv,
            f"exponents_logistic_only_{tag}Hz.png",
            target_freq, source_label,
            logistic_accepted=logistic_accepted
        )
        write_plot_logistic_diag(
            results_dir, legacy_csv,
            f"exponents_logistic_diagnostics_{tag}Hz.png",
            target_freq, source_label,
            logistic_accepted=logistic_accepted
        )

    # Run gnuplot for scripts written in this mode
    scripts = [
        f"rho_vs_S__LOGLOG__{mode}__{methods}_{tag}Hz.gp",
        f"cond_vs_S__LOGLOG__{mode}__{methods}_{tag}Hz.gp",
        f"rho_vs_S__REAL__{mode}__{methods}_{tag}Hz.gp",
        f"cond_vs_S__REAL__{mode}__{methods}_{tag}Hz.gp",
    ]
    for sc in scripts:
        _run_gnuplot(sc, cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)

    if mode == "NORM":
        _run_gnuplot("plot_logistic_only.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)
        _run_gnuplot("plot_logistic_diagnostics.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

    # =========================
    # ALL-VALUES ACCUMULATION (for Excel)
    # =========================
    for _, r in sip_df.iterrows():
        row = dict(r)
        row["run_root"] = str(run_root)
        row["source"] = source_label
        row["mode"] = mode
        row["ici_mode"] = ici_mode
        row["freq_target"] = target_freq
        row["tag"] = tag
        row["logbook_sheet"] = used_sheet
        row["logbook_overlap"] = overlap_cnt
        row["logbook_block_rows"] = block_rows
        all_rows_accum.append(row)
        per_run_excel_rows.append(row)

    return {
        "run": str(run_root),
        "source": source_label,
        "mode": mode,
        "status": "OK",
        "points": int(len(sip_df)),
        "logbook_sheet": used_sheet,
        "logbook_overlap": overlap_cnt,
        "logbook_block_rows": block_rows,
        "poly_deg": int(poly_deg) if (poly_deg and poly_deg >= 2) else 0,
        "log4_enabled": bool(do_log4),
        "logistic_attempted": bool(logistic_attempted),
        "logistic_accepted": bool(logistic_accepted),
        "logistic_status": str(logistic_status),
        "logistic_reason": str(logistic_reason),
        "r2_logrho_lin": float(r2_lr) if np.isfinite(r2_lr) else r2_lr,
        "r2_logsig_lin": float(r2_ls) if np.isfinite(r2_ls) else r2_ls,
        "r2_logrho_poly": float(r2_logrho_poly) if np.isfinite(r2_logrho_poly) else r2_logrho_poly,
        "r2_logsig_poly": float(r2_logsig_poly) if np.isfinite(r2_logsig_poly) else r2_logsig_poly,
        "r2_rho_real_lin": float(r2r) if np.isfinite(r2r) else r2r,
        "r2_sig_real_lin": float(r2s) if np.isfinite(r2s) else r2s,
        "r2_rho_real_poly": float(r2_rho_poly_real) if np.isfinite(r2_rho_poly_real) else r2_rho_poly_real,
        "r2_sig_real_poly": float(r2_sig_poly_real) if np.isfinite(r2_sig_poly_real) else r2_sig_poly_real,
    }

# =============================================================================
# MAIN CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("top", help="Top folder containing run subfolders")
    ap.add_argument("--freq", type=float, default=0.01, help="Target frequency in Hz")
    ap.add_argument("--sample", type=str, default=None, help="Optional substring filter on SIP sheet names")
    ap.add_argument("--min-points", type=int, default=4, help="Minimum matched points to accept a run")
    ap.add_argument("--max-df", type=float, default=0.0, help="Max allowed |f_found - f_target| (0 disables)")
    ap.add_argument("--gnuplot-timeout", type=int, default=30, help="Timeout seconds per gnuplot call")
    ap.add_argument("--clean-results", action="store_true", help="Delete files inside each run's Results/ before writing.")

    ap.add_argument("--ici-mode", type=str, default="sigref_over_sig",
                    choices=["sigref_over_sig","sig_over_sigref"],
                    help="Affects normalized index definition (NORM mode).")

    ap.add_argument("--poly-deg", type=int, default=0, help="Polynomial degree (>=2). 0 disables.")
    ap.add_argument("--poly-raw", action="store_true", help="Kept for compatibility (poly already applied).")

    ap.add_argument("--logistic-fast", action="store_true", help="Use SciPy-free fast logistic search")
    ap.add_argument("--logistic-max-seconds", type=float, default=8.0, help="Time cap for logistic per run")

    ap.add_argument("--logistic4", action="store_true", help="Enable 4-parameter logistic fits (LOG4).")

    ap.add_argument("--both-modes", action="store_true",
                    help="Run BOTH NORM and UNORM modes (recommended).")

    args = ap.parse_args()

    base = Path(args.top).expanduser().resolve()
    if not base.exists():
        print(f"Not found: {base}")
        raise SystemExit(2)

    roots = discover_run_roots_fast(base)
    print(f"Found {len(roots)} run folder(s).", flush=True)

    tag = tag_from_freq(args.freq)

    all_rows = []
    summary_rows = []
    t0 = time.time()

    for i, run in enumerate(roots, 1):
        print(f"[{i}/{len(roots)}] START {run}", flush=True)
        ti = time.time()

        per_run_excel_rows = []

        try:
            modes = ["NORM", "UNORM"] if args.both_modes else ["NORM"]

            for mode in modes:
                res = process_one_mode(
                    run_root=run,
                    target_freq=args.freq,
                    sheet_filter=args.sample,
                    min_points=args.min_points,
                    max_df=args.max_df,
                    do_log4=args.logistic4,
                    logistic_fast=args.logistic_fast,
                    logistic_max_seconds=args.logistic_max_seconds,
                    gnuplot_timeout=args.gnuplot_timeout,
                    clean_results=args.clean_results,
                    ici_mode=args.ici_mode,
                    poly_deg=args.poly_deg,
                    poly_raw=args.poly_raw,
                    mode=mode,
                    all_rows_accum=all_rows,
                    per_run_excel_rows=per_run_excel_rows
                )
                res["i"] = i
                summary_rows.append(res)

            dt = time.time() - ti
            pts = [r.get("points",0) for r in summary_rows if r.get("run")==str(run)]
            pts_txt = "+".join(str(p) for p in pts) if pts else "0"
            mode_txt = "+".join(modes)
            print(f"[{i}/{len(roots)}] OK modes={mode_txt} points={pts_txt} ({dt:.1f}s)", flush=True)

            if per_run_excel_rows:
                df_run = pd.DataFrame(per_run_excel_rows)
                out_xlsx = (run / RESULTS_DIR) / f"all_values_{tag}Hz.xlsx"
                with pd.ExcelWriter(out_xlsx, engine="openpyxl") as w:
                    df_run.to_excel(w, sheet_name="all_values", index=False)

        except subprocess.TimeoutExpired:
            summary_rows.append({"run": str(run), "mode": "?", "status": "FAIL", "reason": "gnuplot timeout", "points": 0})
            print(f"[{i}/{len(roots)}] FAIL gnuplot timeout", flush=True)
        except Exception as e:
            summary_rows.append({"run": str(run), "mode": "?", "status": "FAIL", "reason": str(e), "points": 0})
            print(f"[{i}/{len(roots)}] FAIL {e}", flush=True)

    if summary_rows:
        df_sum = pd.DataFrame(summary_rows)
        sum_csv = base / f"Results_Summary_{tag}Hz.csv"
        sum_xlsx = base / f"Results_Summary_{tag}Hz.xlsx"
        df_sum.to_csv(sum_csv, index=False)
        with pd.ExcelWriter(sum_xlsx, engine="openpyxl") as w:
            df_sum.to_excel(w, sheet_name="summary", index=False)
        print(f"\nSummary written: {sum_csv}", flush=True)
        print(f"Summary written: {sum_xlsx}", flush=True)

    if all_rows:
        df_all = pd.DataFrame(all_rows)
        all_xlsx = base / f"all_values_{tag}Hz.xlsx"
        with pd.ExcelWriter(all_xlsx, engine="openpyxl") as w:
            df_all.to_excel(w, sheet_name="all_values", index=False)
        print(f"All-values written: {all_xlsx}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)

if __name__ == "__main__":
    main()

# Runner you requested:
# python3 auto_sip_pipeline_v27_noscipy.py "/mnt/3rd900/Projects/LA Project new/For_LBNL/SIP" \
#   --freq 0.01 --min-points 4 \
#   --logistic-fast --logistic-max-seconds 8 \
#   --clean-results --both-modes \
#   --ici-mode sig_over_sigref \
#   --poly-deg 2 --poly-raw \
#   --logistic4
