#!/usr/bin/env python3
"""
auto_sip_pipeline_v27_noscipy.py  (UPDATED: FULLY DOCUMENTED / COMMENTED)

v27 + (YOUR NEW REQUESTS ADDED)
===============================

We keep ALL previous approaches and plots EXACTLY as in v27, and we ADD:

A) FULL AUDIT TRAIL TO EXCEL (PER RUN + MASTER)
-----------------------------------------------
For every run folder we write:
  Results/full_audit_<tag>Hz.xlsx
with sheets:
  - data_long   : every raw value, log value, derived calc, fit predictions, residuals
  - fit_params  : all fit parameters + metrics (old + new)
  - notes       : file paths, selected logbook sheet/block, overlap, dropped rows, etc.

At the base folder we write:
  Results_AllRuns_<tag>Hz.xlsx
with sheets:
  - all_data_long  : concatenated "data_long" from all runs (includes run/source columns)
  - all_fit_params : concatenated fit_params (includes run/source columns)

B) REAL-SPACE FITS (NOT LOG PLOTS) + PLOTS
------------------------------------------
We ADD fits in real space (plots in real y vs real x), specifically targeting:
  "electric conductivity (sigma_imag) vs moisture content"
Moisture content:
  - if porosity (phi) is available: theta = S*phi (preferred x)
  - else fallback: x = S

We add NON-NORMALIZED fits (direct values, no ref normalization):
  - sigma_imag vs theta (or S)
  - rho vs theta (or S)
  - optional sigma_real = 1/rho vs theta (or S)

Models (SciPy-free):
  - linear     : y = a*x + c
  - exponential: y = A*exp(b*x)       (requires y>0)
  - power law  : y = A*x^b            (requires x>0 and y>0)

We write all params/preds/residuals to CSV + Excel.

C) MANY NEW PLOTS WITH CLEAR NAMES
----------------------------------
We KEEP old filenames exactly:
  target_exponents_<tag>Hz.png
  index_exponents_<tag>Hz.png
  exponents_raw_<tag>Hz.png
  exponents_logistic_only_<tag>Hz.png
  exponents_logistic_diagnostics_<tag>Hz.png
  plus optional existing v27 poly/vn-smin outputs...

We ADD new plots with explicit searchable naming:
  NONORM_SIGMA_REAL_vs_THETA_ALL_<tag>Hz.png   (preferred if theta exists)
  NONORM_SIGMA_REAL_vs_S_ALL_<tag>Hz.png       (fallback)
  NONORM_RHO_REAL_vs_THETA_ALL_<tag>Hz.png
  NONORM_RHO_REAL_vs_S_ALL_<tag>Hz.png
  (optional) NONORM_SIGMAREAL_REAL_vs_THETA_ALL_<tag>Hz.png  (sigma_real=1/rho)

All new gnuplot scripts use column("name") and titles use the folder source label.
We continue formatting "S-1" as "S_{-1}" via pretty_source_name().

IMPORTANT: "Fits on real values not log"
----------------------------------------
We *plot in real space* (y vs x). For exponential/power models, parameters are estimated
using log-linear regression as a SciPy-free closed-form approach, then evaluated in
real space and reported with REAL-space SSE and REAL-space R^2. This gives you what you
need: real-space curves + real-space metrics + robust SciPy-free fitting.

Final goal
----------
The final goal is to discover a reliable relation between electric conductivity and
moisture content. This update produces exactly the per-run and all-run tables and
plots needed to explore:
  sigma_imag (or sigma_real) vs theta (or S) across all runs.

USAGE EXAMPLES (same as v27, now with extra outputs)
----------------------------------------------------
# baseline
python3 auto_sip_pipeline_v27_noscipy.py "/path/to/SIP" \
  --freq 0.01 --min-points 4 \
  --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
  --clean-results --ici-mode sig_over_sigref

# add VN Smin
python3 auto_sip_pipeline_v27_noscipy.py "/path/to/SIP" \
  --freq 0.01 --min-points 4 \
  --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
  --clean-results --ici-mode sig_over_sigref \
  --vn-smin 0.30

# add polynomial fits on TARGET (deg 2) and RAW logRho
python3 auto_sip_pipeline_v27_noscipy.py "/path/to/SIP" \
  --freq 0.01 --min-points 4 \
  --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
  --clean-results --ici-mode sig_over_sigref \
  --vn-smin 0.30 --poly-deg 2 --poly-raw
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
    """Pick the first matching column name (case-insensitive) from candidates."""
    low = {_norm(c): c for c in df.columns}
    for k in cands:
        kk = _norm(k)
        if kk in low:
            return low[kk]
    return ""

def find_first_matching_file(folder: Path, hint: str, ext=".xlsx"):
    """Find the first file under folder whose name contains 'hint' (case-insensitive) and ends with ext."""
    if not folder.exists():
        return None
    for p in folder.rglob(f"*{ext}"):
        if hint.lower() in p.name.lower():
            return p
    return None

def to_float(x):
    """Robust float conversion (handles '40%' -> 40). Returns NaN on failure."""
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
    """Safe log10: returns NaN if non-finite or <=0."""
    x = float(x)
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return math.log10(x)

def sse(y, yhat):
    """Sum of squared errors."""
    r = y - yhat
    return float(np.sum(r*r))

def _is_blank(x) -> bool:
    """Blank-ish check for logbook cells."""
    if x is None:
        return True
    s = str(x).strip()
    return (s == "") or (s.lower() in ["nan", "none"])

def norm_measid(s: str) -> str:
    """
    Normalize measurement IDs so that 'M 1', 'M-1', 'M_1', 'm1' all match => 'm1'.
    """
    s = str(s).strip().lower()
    s = re.sub(r'[^a-z0-9]+', '', s)
    return s

def pretty_source_name(run_folder_name: str) -> str:
    """
    Plot label from folder name:
      - underscores -> spaces
      - "S-1" -> "S_{-1}" (gnuplot enhanced)
    """
    s = str(run_folder_name).replace("_", " ").strip()
    s = re.sub(r'\bS-(\d+)\b', r'S_{-\1}', s)
    return s

def real_r2(y, yhat):
    """R^2 in real space."""
    y = np.asarray(y, float); yhat = np.asarray(yhat, float)
    m = np.isfinite(y) & np.isfinite(yhat)
    y = y[m]; yhat = yhat[m]
    if len(y) < 2:
        return np.nan
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    return 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan


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
    """Survival function for F distribution."""
    if not np.isfinite(F) or F < 0:
        return np.nan
    x = (d1*F)/(d1*F + d2)
    cdf = betai(d1/2.0, d2/2.0, x)
    return max(0.0, 1.0 - cdf)

def f_test_pvalue(ss_res, ss_tot, df_model, n):
    """F-test p-value for regression models (SciPy-free)."""
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
    """Gaussian AIC."""
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)


# =============================================================================
# FIT ROUTINES: LINEAR / THROUGH ORIGIN / POLYNOMIAL
# =============================================================================

def linear_fit(x, y):
    """
    OLS: y = slope*x + intercept.
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
    """Least squares forced through origin: y = a*x. Returns a, r2, ss_res."""
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
    """Polynomial regression y ~ poly(x,deg). Returns (coeffs, yhat_on_fit_points, r2, ss_res)."""
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
# REAL-SPACE MOISTURE RELATION FITS (NON-NORMALIZED)
# =============================================================================

def fit_linear_real(x, y):
    """Fit y = a*x + c in real space."""
    a, c, r2, _, ss_res, _ = linear_fit(x, y)
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    yhat = (a*x + c) if (np.isfinite(a) and np.isfinite(c)) else None
    return float(a), float(c), yhat, float(r2), float(ss_res)

def fit_exponential_real(x, y):
    """
    Fit y = A * exp(b*x), y>0.
    Parameters estimated via log-linear regression; metrics computed in real space.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y) & (y > 0)
    x = x[m]; y = y[m]
    if len(x) < 3:
        return np.nan, np.nan, None, np.nan, np.nan
    ly = np.log(y)
    A = np.vstack([x, np.ones_like(x)]).T
    b, lnA = np.linalg.lstsq(A, ly, rcond=None)[0]
    A0 = float(np.exp(lnA)); b0 = float(b)
    yhat = A0 * np.exp(b0 * x)
    ss_res = sse(y, yhat)
    r2 = real_r2(y, yhat)
    return A0, b0, yhat, float(r2), float(ss_res)

def fit_powerlaw_real(x, y):
    """
    Fit y = A * x^b, requires x>0 and y>0.
    Parameters estimated via log-log regression; metrics computed in real space.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = x[m]; y = y[m]
    if len(x) < 3:
        return np.nan, np.nan, None, np.nan, np.nan
    lx = np.log(x); ly = np.log(y)
    A = np.vstack([lx, np.ones_like(lx)]).T
    b, lnA = np.linalg.lstsq(A, ly, rcond=None)[0]
    A0 = float(np.exp(lnA)); b0 = float(b)
    yhat = A0 * (x ** b0)
    ss_res = sse(y, yhat)
    r2 = real_r2(y, yhat)
    return A0, b0, yhat, float(r2), float(ss_res)


# =============================================================================
# LOGISTIC FIT (NO SCIPY)
# =============================================================================

def logistic4(x, L, k, x0, c):
    """4-parameter logistic."""
    z = k*(x - x0)
    z = np.clip(z, -60.0, 60.0)
    return c + L/(1.0 + np.exp(-z))

def _logistic_init_guess(x, y):
    """Heuristic init for logistic."""
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
    Fast coordinate grid search around init guess (SciPy-free).
    Returns (fit_dict, status_str).
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    n = len(x)
    if n < 5:
        return None, "not_enough_points"
    if not fast:
        return None, "fail"

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
    """Accept logistic if R2 improves by eps or AIC improves (if enabled)."""
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
      - each sheet = measurement ID
      - take row closest to target_freq
      - extract rho and sigma_imag
    """
    xls = pd.ExcelFile(sip_xlsx)
    out = []
    for sh in xls.sheet_names:
        if sheet_filter and (sheet_filter.lower() not in str(sh).lower()):
            continue

        df = pd.read_excel(sip_xlsx, sheet_name=sh)
        if df.empty:
            continue

        freq_col = pick_col(df, SIP_FREQ_CANDS) or df.columns[0]
        rho_col  = pick_col(df, SIP_RHO_CANDS)
        imag_col = pick_col(df, SIP_IMAG_CANDS)

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
# LOGBOOK READING: ROBUST BLOCK SELECTION (UPDATED: RETURN PHI MAP TOO)
# =============================================================================

def _contiguous_blocks(mask: np.ndarray):
    """Return list of contiguous True blocks (i0,i1) from boolean mask."""
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
    Choose best sheet/block by overlap with SIP measurement IDs.

    Returns:
      sat_map: measurement_norm -> S
      phi_map: measurement_norm -> phi (if available; else missing)
      used_sheet, overlap_cnt, block_rows
    """
    xls = pd.ExcelFile(log_xlsx)
    sip_set = set([m for m in sip_measurements_norm if not _is_blank(m)])

    best = None
    # best tuple:
    # (score, sheet, blk_df, id_col, sat_col, mw_col, por_col, vol_col, overlap)

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
    phi_map = {}

    for _, r in blk.iterrows():
        mid = str(r.get(id_col, "")).strip()
        if _is_blank(mid):
            continue
        midn = norm_measid(mid)

        # Saturation S
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
            sat_map[midn] = float(S)

        # Porosity phi (optional; used for theta = S*phi)
        if por_col:
            phi_val = to_float(r.get(por_col, np.nan))
            if np.isfinite(phi_val) and phi_val > 0:
                phi_map[midn] = float(phi_val)

    if not sat_map:
        raise RuntimeError(f"Selected logbook sheet '{sh}' but mapped 0 saturation rows.")

    return sat_map, phi_map, sh, int(overlap), int(len(blk))


# =============================================================================
# RESULTS CLEANING
# =============================================================================

def clean_results_dir(results_dir: Path):
    """Delete files in Results/ (do not delete folder)."""
    if not results_dir.exists():
        return
    for p in results_dir.iterdir():
        try:
            if p.is_file() or p.is_symlink():
                p.unlink()
        except Exception:
            pass


# =============================================================================
# GNUPLOT SCRIPT WRITERS (OLD v27: KEEP AS-IS)
# =============================================================================
# NOTE: These functions are your existing v27 plot writers and are kept intact.

def write_plot_target(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                      source_label: str,
                      n_vn: float, r2_ri_vn: float, p_vn: float, r2_ici_vn: float):
    label = "{:g}".format(target_freq)
    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set tics out
set key top right
set xlabel 'log10(Saturation)'

set title '{src}: Saturation Exponent (Resistivity Index) at {label} Hz'
set ylabel 'log10(Resistivity Index)'
plot '{csv}' using (column("logS")):(column("logRI_vn")) with points pt 7 ps 1.8 lc rgb '#0066ff' title 'data', \\
     '{csv}' using (column("logS")):(column("logRI_vn_fit")) with lines dt 2 lw 3 lc rgb '#0066ff' \\
     title sprintf('{label}Hz: n = %.2f, R^2 = %.3f', {nval}, {r2a})

set title '{src}: Saturation Exponent (Imag Conductivity Index) at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
plot '{csv}' using (column("logS")):(column("logICI_vn")) with points pt 7 ps 1.8 lc rgb '#dd0000' title 'data', \\
     '{csv}' using (column("logS")):(column("logICI_vn_fit")) with lines dt 2 lw 3 lc rgb '#dd0000' \\
     title sprintf('{label}Hz: p = %.2f, R^2 = %.3f', {pval}, {r2b})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label, src=source_label,
           nval=n_vn, r2a=r2_ri_vn, pval=p_vn, r2b=r2_ici_vn)
    (results_dir / "plot_target_exponents.gp").write_text(gp)

def write_plot_target_smin(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                           source_label: str, vn_smin: float,
                           n_vn_smin: float, r2_ri_vn_smin: float, p_vn_smin: float, r2_ici_vn_smin: float):
    label = "{:g}".format(target_freq)
    if not (np.isfinite(n_vn_smin) and np.isfinite(p_vn_smin) and np.isfinite(r2_ri_vn_smin) and np.isfinite(r2_ici_vn_smin)):
        gp = """reset
set term pngcairo size 1200,700 enhanced font 'Arial,26'
set output '{out_png}'
set grid
unset key
set title '{src}: TARGET (S >= {smin:g}) at {label} Hz'
set label 1 'Not enough points for S >= {smin:g}' at graph 0.5,0.5 center
set xrange [-2:0]
set yrange [-1:1]
plot '-' using 1:2 with points pt 7 ps 0 notitle
-2 0
0  0
e
""".format(out_png=out_png, src=source_label, smin=vn_smin, label=label)
        (results_dir / "plot_target_exponents_smin.gp").write_text(gp)
        return

    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set tics out
set key top right
set xlabel 'log10(Saturation)'

set title '{src}: TARGET (S >= {smin:g}) - Resistivity Index at {label} Hz'
set ylabel 'log10(Resistivity Index)'
plot '{csv}' using ((column("S")>={smin:g})?column("logS"):1/0):((column("S")>={smin:g})?column("logRI_vn"):1/0) \\
     with points pt 7 ps 1.8 lc rgb '#0066ff' title 'data (S>=Smin)', \\
     '{csv}' using ((column("S")>={smin:g})?column("logS"):1/0):((column("S")>={smin:g})?column("logRI_vn_fit_smin"):1/0) \\
     with lines dt 2 lw 3 lc rgb '#0066ff' \\
     title sprintf('{label}Hz: n = %.2f, R^2 = %.3f', {nval}, {r2a})

set title '{src}: TARGET (S >= {smin:g}) - Imag Conductivity Index at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
plot '{csv}' using ((column("S")>={smin:g})?column("logS"):1/0):((column("S")>={smin:g})?column("logICI_vn"):1/0) \\
     with points pt 7 ps 1.8 lc rgb '#dd0000' title 'data (S>=Smin)', \\
     '{csv}' using ((column("S")>={smin:g})?column("logS"):1/0):((column("S")>={smin:g})?column("logICI_vn_fit_smin"):1/0) \\
     with lines dt 2 lw 3 lc rgb '#dd0000' \\
     title sprintf('{label}Hz: p = %.2f, R^2 = %.3f', {pval}, {r2b})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label, src=source_label, smin=vn_smin,
           nval=n_vn_smin, r2a=r2_ri_vn_smin, pval=p_vn_smin, r2b=r2_ici_vn_smin)

    (results_dir / "plot_target_exponents_smin.gp").write_text(gp)

def write_plot_target_poly(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                           source_label: str, deg: int,
                           r2_ri_poly: float, r2_ici_poly: float):
    label = "{:g}".format(target_freq)
    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set tics out
set key top right
set xlabel 'log10(Saturation)'

set title '{src}: TARGET polynomial (deg {deg}) - Resistivity Index at {label} Hz'
set ylabel 'log10(Resistivity Index)'
plot '{csv}' using (column("logS")):(column("logRI_vn")) with points pt 7 ps 1.8 lc rgb '#0066ff' title 'data', \\
     '{csv}' using (column("logS")):(column("logRI_vn_polyfit")) with lines dt 1 lw 4 lc rgb '#0066ff' \\
     title sprintf('poly deg {deg}: R^2 = %.3f', {r2a})

set title '{src}: TARGET polynomial (deg {deg}) - Imag Conductivity Index at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
plot '{csv}' using (column("logS")):(column("logICI_vn")) with points pt 7 ps 1.8 lc rgb '#dd0000' title 'data', \\
     '{csv}' using (column("logS")):(column("logICI_vn_polyfit")) with lines dt 1 lw 4 lc rgb '#dd0000' \\
     title sprintf('poly deg {deg}: R^2 = %.3f', {r2b})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label, src=source_label, deg=deg,
           r2a=r2_ri_poly if np.isfinite(r2_ri_poly) else float("nan"),
           r2b=r2_ici_poly if np.isfinite(r2_ici_poly) else float("nan"))
    (results_dir / "plot_target_exponents_poly.gp").write_text(gp)

def write_plot_target_poly_smin(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                                source_label: str, deg: int, vn_smin: float,
                                r2_ri_poly: float, r2_ici_poly: float):
    label = "{:g}".format(target_freq)
    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set tics out
set key top right
set xlabel 'log10(Saturation)'

set title '{src}: TARGET poly (deg {deg}, S >= {smin:g}) - Resistivity Index at {label} Hz'
set ylabel 'log10(Resistivity Index)'
plot '{csv}' using ((column("S")>={smin:g})?column("logS"):1/0):((column("S")>={smin:g})?column("logRI_vn"):1/0) \\
     with points pt 7 ps 1.8 lc rgb '#0066ff' title 'data (S>=Smin)', \\
     '{csv}' using ((column("S")>={smin:g})?column("logS"):1/0):((column("S")>={smin:g})?column("logRI_vn_polyfit_smin"):1/0) \\
     with lines dt 1 lw 4 lc rgb '#0066ff' \\
     title sprintf('poly deg {deg}: R^2 = %.3f', {r2a})

set title '{src}: TARGET poly (deg {deg}, S >= {smin:g}) - Imag Conductivity Index at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
plot '{csv}' using ((column("S")>={smin:g})?column("logS"):1/0):((column("S")>={smin:g})?column("logICI_vn"):1/0) \\
     with points pt 7 ps 1.8 lc rgb '#dd0000' title 'data (S>=Smin)', \\
     '{csv}' using ((column("S")>={smin:g})?column("logS"):1/0):((column("S")>={smin:g})?column("logICI_vn_polyfit_smin"):1/0) \\
     with lines dt 1 lw 4 lc rgb '#dd0000' \\
     title sprintf('poly deg {deg}: R^2 = %.3f', {r2b})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label, src=source_label, deg=deg, smin=vn_smin,
           r2a=r2_ri_poly if np.isfinite(r2_ri_poly) else float("nan"),
           r2b=r2_ici_poly if np.isfinite(r2_ici_poly) else float("nan"))
    (results_dir / "plot_target_exponents_poly_smin.gp").write_text(gp)

def write_plot_index(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                     source_label: str,
                     n_idx: float, r2_ri: float, p_idx: float, r2_ici: float):
    label = "{:g}".format(target_freq)
    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set tics out
set key top right
set xlabel 'log10(Saturation)'

set title '{src}: Resistivity Index vs Saturation (goal, through-origin) at {label} Hz'
set ylabel 'log10(Resistivity Index)'
a1 = {a1}
f1(x) = a1*x
plot '{csv}' using (column("logS")):(column("logRhoIdx")) with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('{label}Hz: n=%.2f, R^2=%.3f', {n_idx}, {r2_ri})

set title '{src}: Imag Conductivity Index vs Saturation (goal, through-origin) at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
a2 = {a2}
f2(x) = a2*x
plot '{csv}' using (column("logS")):(column("logImagIdx")) with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('{label}Hz: p=%.2f, R^2=%.3f', {p_idx}, {r2_ici})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label, src=source_label,
           a1=(-n_idx), n_idx=n_idx, r2_ri=r2_ri,
           a2=(-p_idx), p_idx=p_idx, r2_ici=r2_ici)
    (results_dir / "plot_index_exponents.gp").write_text(gp)

def write_plot_raw(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                   source_label: str,
                   lin_rho, lin_im, log_fit_acc):
    label = "{:g}".format(target_freq)
    a_r = lin_rho["slope"]; b_r = lin_rho["intercept"]; r2_r = lin_rho["r2"]; p_r = lin_rho["p"]
    a_i = lin_im["slope"];  b_i = lin_im["intercept"];  r2_i = lin_im["r2"];  p_i = lin_im["p"]

    log_defs = ""
    log_clause = ""
    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        r2_log = log_fit_acc.get("r2", np.nan)
        p_model = log_fit_acc.get("p_model", np.nan)
        log_defs = (
            "L1={:.12g}\n"
            "k1={:.12g}\n"
            "x01={:.12g}\n"
            "c1={:.12g}\n"
            "fL(x)=c1+L1/(1+exp(-k1*(x-x01)))\n"
        ).format(L, k, x0, c)
        if np.isfinite(r2_log) and np.isfinite(p_model):
            log_clause = ", \\\n     fL(x) with lines dt 1 lw 4 title sprintf('logistic: R^2=%.3f, p=%.3g', {:.12g}, {:.12g})".format(r2_log, p_model)
        else:
            log_clause = ", \\\n     fL(x) with lines dt 1 lw 4 title 'logistic fit'"

    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set key top right

set title '{src}: Resistivity vs Saturation (log-log) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
a1={a1:.12g}
b1={b1:.12g}
f1(x)=a1*x+b1
{log_defs}
plot '{csv}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('linear: n=%.2f, R^2=%.3f, p=%.3g', -a1, {r2_r:.12g}, {p_r:.12g}){log_clause}

set title '{src}: Imag Conductivity vs Saturation (log-log) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Imag Conductivity)'
a2={a2:.12g}
b2={b2:.12g}
f2(x)=a2*x+b2
plot '{csv}' using (column("logS")):(column("logImag")) with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('linear: p=%.2f, R^2=%.3f, p=%.3g', a2, {r2_i:.12g}, {p_i:.12g})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label, src=source_label,
           a1=a_r, b1=b_r, r2_r=r2_r, p_r=p_r,
           a2=a_i, b2=b_i, r2_i=r2_i, p_i=p_i,
           log_defs=log_defs, log_clause=log_clause)
    (results_dir / "plot_raw_exponents.gp").write_text(gp)

def write_plot_raw_poly(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                        source_label: str, deg: int, r2_poly: float):
    label = "{:g}".format(target_freq)
    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,26'
set output '{out_png}'

set grid
set key top right
set title '{src}: RAW polynomial (deg {deg}) - logRho vs logS at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'

plot '{csv}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.7 title 'data', \\
     '{csv}' using (column("logS")):(column("logRho_polyfit")) with lines dt 1 lw 4 title sprintf('poly deg {deg}: R^2 = %.3f', {r2})
""".format(out_png=out_png, csv=csv_name, src=source_label, label=label, deg=deg,
           r2=r2_poly if np.isfinite(r2_poly) else float("nan"))
    (results_dir / "plot_raw_exponents_poly.gp").write_text(gp)

def write_plot_logistic_only(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                            source_label: str,
                            log_fit_acc):
    label = "{:g}".format(target_freq)
    log_defs = ""
    log_line = ""
    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        r2_log = log_fit_acc.get("r2", np.nan)
        p_model = log_fit_acc.get("p_model", np.nan)
        log_defs = (
            "L1={:.12g}\n"
            "k1={:.12g}\n"
            "x01={:.12g}\n"
            "c1={:.12g}\n"
            "fL(x)=c1+L1/(1+exp(-k1*(x-x01)))\n"
        ).format(L, k, x0, c)
        if np.isfinite(r2_log) and np.isfinite(p_model):
            log_line = "fL(x) with lines dt 1 lw 4 title sprintf('logistic: R^2=%.3f, p=%.3g', {:.12g}, {:.12g})".format(r2_log, p_model)
        else:
            log_line = "fL(x) with lines dt 1 lw 4 title 'logistic fit'"

    if not log_line:
        log_line = "1/0 notitle"

    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,26'
set output '{out_png}'

set title '{src}: Resistivity vs Saturation (log-log) at {label} Hz (logistic only)'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'

{log_defs}
plot '{csv}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.7 title 'data', \\
     {log_line}
""".format(out_png=out_png, csv=csv_name, label=label, src=source_label, log_defs=log_defs, log_line=log_line)
    (results_dir / "plot_logistic_only.gp").write_text(gp)

def write_plot_logistic_diag(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                             source_label: str,
                             has_log: bool):
    label = "{:g}".format(target_freq)
    extra_curve = ""
    if has_log:
        extra_curve = ", \\\n     '{csv}' using (column(\"logS\")):(column(\"logRho_logfit\")) with lines dt 1 lw 4 title 'logistic (accepted)'".format(csv=csv_name)

    if has_log:
        panel3 = "plot '{csv}' using (column(\"logS\")):(column(\"resid_log_rho\")) with points pt 7 ps 1.4 title 'residuals'\n".format(csv=csv_name)
        extra3 = ""
    else:
        extra3 = (
            "unset key\n"
            "set label 1 'No logistic accepted' at graph 0.5,0.5 center front\n"
            "set xrange [-2:0]\n"
            "set yrange [-1:1]\n"
        )
        panel3 = (
            "plot '-' using 1:2 with points pt 7 ps 0 notitle\n"
            "-2 0\n"
            "0  0\n"
            "e\n"
        )

    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1800 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 3,1

set title '{src}: Diagnostics: logRho vs logS at {label} Hz'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
plot '{csv}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.6 title 'data', \\
     '{csv}' using (column("logS")):(column("logRho_linfit")) with lines dt 2 lw 3 title 'linear'{extra_curve}

set title 'Residuals: linear (data - linear)'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'residual'
plot '{csv}' using (column("logS")):(column("resid_lin_rho")) with points pt 7 ps 1.4 title 'residuals'

set title 'Residuals: logistic (data - logistic)'
set grid
set xlabel 'log10(Saturation)'
set ylabel 'residual'
{extra3}
{panel3}

unset multiplot
""".format(out_png=out_png, label=label, csv=csv_name, src=source_label,
           extra_curve=extra_curve, extra3=extra3, panel3=panel3)
    (results_dir / "plot_logistic_diagnostics.gp").write_text(gp)


# =============================================================================
# GNUPLOT WRITER FOR NEW NON-NORMALIZED REAL-SPACE PLOTS (NEW)
# =============================================================================

def write_plot_nonorm_real(results_dir: Path,
                           csv_name: str,
                           out_png: str,
                           target_freq: float,
                           source_label: str,
                           x_col: str,
                           y_col: str,
                           y_fit_cols: dict,
                           r2_map: dict):
    """
    Plot real-space (non-normalized) relation y vs x.

    x_col: e.g. "theta" or "S"
    y_col: e.g. "sigma_imag" or "rho" or "sigma_real"

    y_fit_cols: dict name->csv_column, for overlays (e.g. {"LIN":"sigma_nonorm_fit_lin", ...})
    r2_map: dict name->r2 value for title strings

    The plot overlays all available fits; missing fit cols are plotted as 1/0.
    """
    label = "{:g}".format(target_freq)
    xlab = x_col
    ylab = y_col
    title = f"{source_label}: NONORM {y_col} vs {x_col} (real space) at {label} Hz"

    # Build plot clauses for fits
    fit_clauses = []
    for name, col in y_fit_cols.items():
        r2v = r2_map.get(name, np.nan)
        if np.isfinite(r2v):
            fit_clauses.append(
                "'{csv}' using (column(\"{x}\")):(column(\"{y}\")) with lines lw 4 title sprintf('{nm} R^2=%.3f', {r2})"
                .format(csv=csv_name, x=x_col, y=col, nm=name, r2=float(r2v))
            )
        else:
            fit_clauses.append(
                "'{csv}' using (column(\"{x}\")):(column(\"{y}\")) with lines lw 4 title '{nm}'"
                .format(csv=csv_name, x=x_col, y=col, nm=name)
            )

    fit_block = ""
    if fit_clauses:
        fit_block = ", \\\n     " + ", \\\n     ".join(fit_clauses)

    gp = """reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1400,900 enhanced font 'Arial,26'
set output '{out_png}'

set grid
set tics out
set key top left
set title '{title}'
set xlabel '{xlab}'
set ylabel '{ylab}'

plot '{csv}' using (column("{xcol}")):(column("{ycol}")) with points pt 7 ps 1.8 title 'data'{fit_block}
""".format(out_png=out_png, title=title, xlab=xlab, ylab=ylab, csv=csv_name,
           xcol=x_col, ycol=y_col, fit_block=fit_block)

    (results_dir / f"plot_{Path(out_png).stem}.gp").write_text(gp)


# =============================================================================
# DISCOVERY + GNUPLOT EXECUTION
# =============================================================================

def discover_run_roots_fast(base: Path):
    """Run root = folder that contains Analysis/ and SIP_Raw_Data/."""
    roots = set()
    for a in base.rglob(ANALYSIS_DIR):
        if not a.is_dir():
            continue
        run = a.parent
        if (run / RAW_DIR).is_dir():
            roots.add(run)
    return sorted(roots)

def _run_gnuplot(script_name: str, cwd: Path, timeout_s: int, fatal: bool):
    """Execute gnuplot script in cwd."""
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
# MAIN PROCESSING FOR ONE RUN (UPDATED)
# =============================================================================

def process_run(run_root: Path,
                target_freq: float,
                sheet_filter: str|None,
                dry_run: bool,
                min_points: int,
                max_df: float,
                do_logistic: bool,
                logistic_min_points: int,
                logistic_fast: bool,
                logistic_max_seconds: float,
                gnuplot_timeout: int,
                clean_results: bool,
                ici_mode: str,
                vn_smin: float,
                poly_deg: int,
                poly_raw: bool):
    """
    v27 processing + NEW:
      - compute phi + theta if possible
      - compute NONORM real-space fits for sigma/rho/sigma_real vs theta/S
      - write per-run Excel full audit
      - write new NONORM plots with explicit names
    """
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return None

    source_label = pretty_source_name(run_root.name)

    results_dir.mkdir(parents=True, exist_ok=True)
    if clean_results:
        clean_results_dir(results_dir)

    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df, sheet_filter=sheet_filter)
    tag = str(target_freq).replace(".","p")

    if sip_df.empty:
        return {"run": str(run_root), "status": "FAIL",
                "reason": "No SIP sheets readable / target freq missing",
                "points": 0}

    # UPDATED: sat_map + phi_map
    sat_map, phi_map, used_sheet, overlap_cnt, block_rows = read_logbook_bestblock_by_overlap(
        log_xlsx, sip_df["measurement_norm"].tolist()
    )

    sip_df["S"] = sip_df["measurement_norm"].map(sat_map).astype(float)

    # NEW: map phi and compute theta
    sip_df["phi"] = sip_df["measurement_norm"].map(phi_map).astype(float)
    sip_df["theta"] = sip_df["S"] * sip_df["phi"]

    # Compute logs for legacy v27 methods
    sip_df["logS"]    = sip_df["S"].apply(safe_log10)
    sip_df["logRho"]  = sip_df["rho"].apply(safe_log10)
    sip_df["logImag"] = sip_df["sigma_imag"].apply(safe_log10)

    # Optional "sigma_real" from resistivity (non-normalized derived conductivity)
    sip_df["sigma_real"] = np.where((np.isfinite(sip_df["rho"])) & (sip_df["rho"] > 0), 1.0/sip_df["rho"], np.nan)

    dbg_all = sip_df.copy()
    dbg_all.to_csv(results_dir / f"debug_before_filter_{tag}Hz.csv", index=False)

    mask = np.isfinite(sip_df["logS"]) & np.isfinite(sip_df["logRho"]) & np.isfinite(sip_df["logImag"])
    dropped = sip_df.loc[~mask, ["measurement", "S", "rho", "sigma_imag", "phi", "theta"]].copy()
    if len(dropped) > 0:
        dropped.to_csv(results_dir / f"dropped_rows_{tag}Hz.csv", index=False)
    sip_df = sip_df.loc[mask].copy()

    if len(sip_df) < min_points:
        return {"run": str(run_root), "status": "FAIL",
                "reason": f"Not enough matched points (need >= {min_points})",
                "points": int(len(sip_df)),
                "logbook_sheet": used_sheet,
                "logbook_overlap": overlap_cnt,
                "logbook_block_rows": block_rows,
                "source": source_label}

    # Reference point = max S (legacy)
    iref = int(np.nanargmax(sip_df["S"].to_numpy(float)))
    rho_ref  = float(sip_df.iloc[iref]["rho"])
    sig_ref  = float(sip_df.iloc[iref]["sigma_imag"])
    logRho_ref = safe_log10(rho_ref)
    logSig_ref = safe_log10(sig_ref)

    # GOAL/INDEX (legacy)
    sip_df["rho_idx"]   = sip_df["rho"] / rho_ref
    sip_df["logRhoIdx"] = sip_df["logRho"] - logRho_ref

    if ici_mode == "sig_over_sigref":
        sip_df["imag_idx"]   = sip_df["sigma_imag"] / sig_ref
        sip_df["logImagIdx"] = sip_df["logImag"] - logSig_ref
    else:
        sip_df["imag_idx"]   = sig_ref / sip_df["sigma_imag"]
        sip_df["logImagIdx"] = logSig_ref - sip_df["logImag"]

    # RAW fits (legacy)
    slope_r, intercept_r, r2_r, p_r, ss_res_lin_r, aic_lin_r = linear_fit(sip_df["logS"], sip_df["logRho"])
    slope_i, intercept_i, r2_i, p_i, ss_res_lin_i, aic_lin_i = linear_fit(sip_df["logS"], sip_df["logImag"])
    lin_rho = {"slope": slope_r, "intercept": intercept_r, "r2": r2_r, "p": p_r, "aic": aic_lin_r}
    lin_im  = {"slope": slope_i, "intercept": intercept_i, "r2": r2_i, "p": p_i, "aic": aic_lin_i}

    # GOAL/INDEX through-origin (legacy)
    a_ri,  r2_ri,  _ = linear_fit_through_origin(sip_df["logS"], sip_df["logRhoIdx"])
    a_ici, r2_ici, _ = linear_fit_through_origin(sip_df["logS"], sip_df["logImagIdx"])

    n_idx = -a_ri
    if ici_mode == "sig_over_sigref":
        p_idx = +a_ici
    else:
        p_idx = -a_ici

    # TARGET VN style (legacy: always downtrend, intercept)
    sip_df["logRI_vn"]  = sip_df["logRhoIdx"]
    sip_df["logICI_vn"] = logSig_ref - sip_df["logImag"]

    slope_RI_vn,  int_RI_vn,  r2_RI_vn,  _, _, _ = linear_fit(sip_df["logS"], sip_df["logRI_vn"])
    slope_ICI_vn, int_ICI_vn, r2_ICI_vn, _, _, _ = linear_fit(sip_df["logS"], sip_df["logICI_vn"])

    n_vn = -slope_RI_vn
    p_vn = -slope_ICI_vn

    sip_df["logRI_vn_fit"]  = slope_RI_vn  * sip_df["logS"] + int_RI_vn
    sip_df["logICI_vn_fit"] = slope_ICI_vn * sip_df["logS"] + int_ICI_vn

    # VN Smin subset (legacy optional)
    n_vn_smin = np.nan
    p_vn_smin = np.nan
    r2_RI_vn_smin = np.nan
    r2_ICI_vn_smin = np.nan
    sip_df["logRI_vn_fit_smin"]  = np.nan
    sip_df["logICI_vn_fit_smin"] = np.nan

    if vn_smin and vn_smin > 0:
        sub = sip_df[sip_df["S"] >= vn_smin].copy()
        if len(sub) >= min_points:
            slope_RI_vn_smin, int_RI_vn_smin, r2_RI_vn_smin, _, _, _ = linear_fit(sub["logS"], sub["logRI_vn"])
            slope_ICI_vn_smin, int_ICI_vn_smin, r2_ICI_vn_smin, _, _, _ = linear_fit(sub["logS"], sub["logICI_vn"])
            n_vn_smin = -slope_RI_vn_smin
            p_vn_smin = -slope_ICI_vn_smin
            m = sip_df["S"] >= vn_smin
            sip_df.loc[m, "logRI_vn_fit_smin"]  = slope_RI_vn_smin  * sip_df.loc[m, "logS"] + int_RI_vn_smin
            sip_df.loc[m, "logICI_vn_fit_smin"] = slope_ICI_vn_smin * sip_df.loc[m, "logS"] + int_ICI_vn_smin

    # Polynomial TARGET (legacy optional)
    r2_RI_vn_poly = np.nan
    r2_ICI_vn_poly = np.nan
    sip_df["logRI_vn_polyfit"]  = np.nan
    sip_df["logICI_vn_polyfit"] = np.nan

    r2_RI_vn_poly_smin = np.nan
    r2_ICI_vn_poly_smin = np.nan
    sip_df["logRI_vn_polyfit_smin"]  = np.nan
    sip_df["logICI_vn_polyfit_smin"] = np.nan

    if poly_deg and poly_deg >= 2:
        c1, _, r2a, _ = poly_fit_predict(sip_df["logS"], sip_df["logRI_vn"], poly_deg)
        c2, _, r2b, _ = poly_fit_predict(sip_df["logS"], sip_df["logICI_vn"], poly_deg)
        r2_RI_vn_poly = r2a
        r2_ICI_vn_poly = r2b
        if c1 is not None:
            sip_df["logRI_vn_polyfit"] = np.polyval(c1, sip_df["logS"].to_numpy(float))
        if c2 is not None:
            sip_df["logICI_vn_polyfit"] = np.polyval(c2, sip_df["logS"].to_numpy(float))

        if vn_smin and vn_smin > 0:
            sub = sip_df[sip_df["S"] >= vn_smin].copy()
            if len(sub) >= max(min_points, poly_deg + 1):
                c1s, _, r2as, _ = poly_fit_predict(sub["logS"], sub["logRI_vn"], poly_deg)
                c2s, _, r2bs, _ = poly_fit_predict(sub["logS"], sub["logICI_vn"], poly_deg)
                r2_RI_vn_poly_smin = r2as
                r2_ICI_vn_poly_smin = r2bs
                m = sip_df["S"] >= vn_smin
                if c1s is not None:
                    sip_df.loc[m, "logRI_vn_polyfit_smin"] = np.polyval(c1s, sip_df.loc[m, "logS"].to_numpy(float))
                if c2s is not None:
                    sip_df.loc[m, "logICI_vn_polyfit_smin"] = np.polyval(c2s, sip_df.loc[m, "logS"].to_numpy(float))

    # Polynomial RAW logRho (legacy optional)
    r2_raw_poly = np.nan
    sip_df["logRho_polyfit"] = np.nan
    if poly_raw and (poly_deg and poly_deg >= 2):
        cR, _, r2p, _ = poly_fit_predict(sip_df["logS"], sip_df["logRho"], poly_deg)
        r2_raw_poly = r2p
        if cR is not None:
            sip_df["logRho_polyfit"] = np.polyval(cR, sip_df["logS"].to_numpy(float))

    # Logistic on RAW logRho (legacy optional)
    log_fit = None
    log_status = "skipped"
    if do_logistic and len(sip_df) >= logistic_min_points:
        log_fit, log_status = logistic_fit_noscipy(
            sip_df["logS"], sip_df["logRho"],
            fast=logistic_fast,
            max_seconds=(logistic_max_seconds if logistic_max_seconds and logistic_max_seconds > 0 else 0.0)
        )
    elif do_logistic:
        log_status = "not_enough_points"

    log_fit_acc = None
    log_accept_reason = "n/a"
    if log_fit is not None:
        ok, reason = accept_logistic(r2_r, aic_lin_r, log_fit)
        log_accept_reason = reason
        if ok:
            log_fit_acc = log_fit

    # diagnostics residuals (legacy)
    logS_arr = sip_df["logS"].to_numpy(float)
    logR_arr = sip_df["logRho"].to_numpy(float)
    lin_pred = slope_r*logS_arr + intercept_r
    sip_df["resid_lin_rho"] = logR_arr - lin_pred

    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        log_pred = logistic4(logS_arr, L, k, x0, c)
        sip_df["resid_log_rho"] = logR_arr - log_pred
        sip_df["logRho_linfit"] = lin_pred
        sip_df["logRho_logfit"] = log_pred
    else:
        sip_df["resid_log_rho"] = np.nan
        sip_df["logRho_linfit"] = lin_pred
        sip_df["logRho_logfit"] = np.nan

    # -------------------------------------------------------------------------
    # NEW: NON-NORMALIZED REAL-SPACE MOISTURE RELATION FITS
    # -------------------------------------------------------------------------
    # Prefer theta if it has enough finite values; else fallback to S.
    use_theta = np.isfinite(sip_df["theta"]).sum() >= min_points
    x_col = "theta" if use_theta else "S"
    x = sip_df[x_col].to_numpy(float)

    # Helper to fit y vs x and write predictions aligned to the dataframe
    def _apply_nonorm_fits(y_col_name: str, prefix: str):
        y = sip_df[y_col_name].to_numpy(float)

        # Prepare output columns
        col_lin = f"{prefix}_fit_lin"
        col_exp = f"{prefix}_fit_exp"
        col_pow = f"{prefix}_fit_pow"
        sip_df[col_lin] = np.nan
        sip_df[col_exp] = np.nan
        sip_df[col_pow] = np.nan
        sip_df[f"resid_{prefix}_lin"] = np.nan
        sip_df[f"resid_{prefix}_exp"] = np.nan
        sip_df[f"resid_{prefix}_pow"] = np.nan

        # Fit only where finite
        m = np.isfinite(x) & np.isfinite(y)
        xx = x[m]; yy = y[m]
        if len(xx) < min_points:
            return {
                "x_col": x_col, "y_col": y_col_name,
                "lin": None, "exp": None, "pow": None
            }

        # Linear
        a, c0, yhat_lin, r2_lin, sse_lin = fit_linear_real(xx, yy)
        if yhat_lin is not None:
            sip_df.loc[m, col_lin] = yhat_lin
            sip_df.loc[m, f"resid_{prefix}_lin"] = yy - yhat_lin

        # Exponential
        Aexp, bexp, yhat_exp, r2_exp, sse_exp = fit_exponential_real(xx, yy)
        if yhat_exp is not None:
            sip_df.loc[m, col_exp] = yhat_exp
            sip_df.loc[m, f"resid_{prefix}_exp"] = yy - yhat_exp

        # Power law
        Apow, bpow, yhat_pow, r2_pow, sse_pow = fit_powerlaw_real(xx, yy)
        if yhat_pow is not None:
            sip_df.loc[m, col_pow] = yhat_pow
            sip_df.loc[m, f"resid_{prefix}_pow"] = yy - yhat_pow

        return {
            "x_col": x_col, "y_col": y_col_name,
            "lin": {"a": a, "c": c0, "r2": r2_lin, "sse": sse_lin},
            "exp": {"A": Aexp, "b": bexp, "r2": r2_exp, "sse": sse_exp},
            "pow": {"A": Apow, "b": bpow, "r2": r2_pow, "sse": sse_pow},
        }

    # Apply to sigma_imag, rho, and sigma_real
    fit_nonorm_sigma = _apply_nonorm_fits("sigma_imag", "sigma_nonorm")
    fit_nonorm_rho   = _apply_nonorm_fits("rho",       "rho_nonorm")
    fit_nonorm_sigR  = _apply_nonorm_fits("sigma_real","sigma_real_nonorm")

    # -------------------------------------------------------------------------
    # WRITE joined CSV (UPDATED: includes everything old + new nonorm fields)
    # -------------------------------------------------------------------------
    csv_path = results_dir / f"joined_{tag}Hz.csv"

    # Ensure all columns exist even if some fits not possible
    want_cols = [
        "measurement","measurement_norm","freq_hz",
        "S","phi","theta",
        "rho","sigma_imag","sigma_real",
        "logS","logRho","logImag",
        "rho_idx","imag_idx","logRhoIdx","logImagIdx",

        # TARGET VN
        "logRI_vn","logICI_vn","logRI_vn_fit","logICI_vn_fit",
        "logRI_vn_fit_smin","logICI_vn_fit_smin",

        # TARGET poly
        "logRI_vn_polyfit","logICI_vn_polyfit",
        "logRI_vn_polyfit_smin","logICI_vn_polyfit_smin",

        # RAW poly
        "logRho_polyfit",

        # diagnostics
        "resid_lin_rho","resid_log_rho",
        "logRho_linfit","logRho_logfit",

        # NEW nonorm predictions/residuals
        "sigma_nonorm_fit_lin","sigma_nonorm_fit_exp","sigma_nonorm_fit_pow",
        "resid_sigma_nonorm_lin","resid_sigma_nonorm_exp","resid_sigma_nonorm_pow",

        "rho_nonorm_fit_lin","rho_nonorm_fit_exp","rho_nonorm_fit_pow",
        "resid_rho_nonorm_lin","resid_rho_nonorm_exp","resid_rho_nonorm_pow",

        "sigma_real_nonorm_fit_lin","sigma_real_nonorm_fit_exp","sigma_real_nonorm_fit_pow",
        "resid_sigma_real_nonorm_lin","resid_sigma_real_nonorm_exp","resid_sigma_real_nonorm_pow",
    ]

    for c in want_cols:
        if c not in sip_df.columns:
            sip_df[c] = np.nan

    sip_df[want_cols].to_csv(csv_path, index=False)

    # -------------------------------------------------------------------------
    # REPORT (UPDATED: include nonorm fits summary)
    # -------------------------------------------------------------------------
    report = results_dir / f"fit_report_{tag}Hz.txt"
    base_text = (
        f"Run root: {run_root}\n"
        f"Source (folder): {source_label}\n"
        f"SIP file: {sip_xlsx}\n"
        f"Logbook file: {log_xlsx}\n"
        f"Logbook sheet used: {used_sheet}\n"
        f"Logbook block rows (usable): {block_rows}\n"
        f"Overlap with SIP measurement IDs: {overlap_cnt}\n\n"
        f"Matched points used: {len(sip_df)}\n"
        f"Reference point (max S): {sip_df.iloc[iref]['measurement']}  S_ref={sip_df.iloc[iref]['S']:.6g}\n"
        f"rho_ref={rho_ref:.6g}, sigma_ref={sig_ref:.6g}\n\n"
        f"GOAL/INDEX ici_mode: {ici_mode}\n"
        f"GOAL defs: RI=rho/rho_ref\n"
        f"          ICI={'sigma_imag/sigma_ref' if ici_mode=='sig_over_sigref' else 'sigma_ref/sigma_imag'}\n\n"
        f"GOAL (through-origin) exponents:\n"
        f"  n(index)={n_idx:.6g} (R2={r2_ri:.6g})\n"
        f"  p(index)={p_idx:.6g} (R2={r2_ici:.6g})\n\n"
        f"TARGET (Van Nuys style; WITH intercept; both downtrend):\n"
        f"  n_vn={n_vn:.6g} (R2={r2_RI_vn:.6g})\n"
        f"  p_vn={p_vn:.6g} (R2={r2_ICI_vn:.6g})\n\n"
        f"RAW linear:\n"
        f"  logRho=a*logS+b: a={slope_r:.6g}, b={intercept_r:.6g}, n(raw)={-slope_r:.6g}, R2={r2_r:.6g}\n"
        f"  logImag=p*logS+q: p={slope_i:.6g}, q={intercept_i:.6g}, R2={r2_i:.6g}\n\n"
        f"Logistic: status={log_status}, accepted={bool(log_fit_acc is not None)}, reason={log_accept_reason}\n\n"
        f"NONORM REAL-SPACE fits (conductivity/moisture goal):\n"
        f"  x used = {x_col} (theta preferred if available; else S)\n"
    )
    report.write_text(base_text)

    def _append_nonorm_summary(label_y: str, fitdict):
        lines = [f"\n  {label_y} vs {fitdict.get('x_col')}:\n"]
        for nm in ["lin","exp","pow"]:
            f0 = fitdict.get(nm, None)
            if f0 is None:
                lines.append(f"    {nm}: not fit (not enough points or invalid domain)\n")
            else:
                if nm == "lin":
                    lines.append(f"    lin: a={f0['a']:.6g}, c={f0['c']:.6g}, R2_real={f0['r2']:.6g}, SSE_real={f0['sse']:.6g}\n")
                else:
                    lines.append(f"    {nm}: A={f0['A']:.6g}, b={f0['b']:.6g}, R2_real={f0['r2']:.6g}, SSE_real={f0['sse']:.6g}\n")
        report.write_text(report.read_text() + "".join(lines))

    _append_nonorm_summary("sigma_imag", fit_nonorm_sigma)
    _append_nonorm_summary("rho", fit_nonorm_rho)
    _append_nonorm_summary("sigma_real=1/rho", fit_nonorm_sigR)

    if vn_smin and vn_smin > 0:
        report.write_text(report.read_text() +
            f"\nEXTRA TARGET (VN style; WITH intercept; S >= {vn_smin:g}):\n"
            f"  n_vn_smin={n_vn_smin:.6g} (R2={r2_RI_vn_smin:.6g})\n"
            f"  p_vn_smin={p_vn_smin:.6g} (R2={r2_ICI_vn_smin:.6g})\n"
        )

    if poly_deg and poly_deg >= 2:
        report.write_text(report.read_text() +
            f"\nPOLY TARGET (deg {poly_deg}):\n"
            f"  R2_RI_vn_poly={r2_RI_vn_poly:.6g}\n"
            f"  R2_ICI_vn_poly={r2_ICI_vn_poly:.6g}\n"
        )
        if vn_smin and vn_smin > 0:
            report.write_text(report.read_text() +
                f"POLY TARGET (deg {poly_deg}, S >= {vn_smin:g}):\n"
                f"  R2_RI_vn_poly_smin={r2_RI_vn_poly_smin:.6g}\n"
                f"  R2_ICI_vn_poly_smin={r2_ICI_vn_poly_smin:.6g}\n"
            )
        if poly_raw:
            report.write_text(report.read_text() +
                f"POLY RAW logRho (deg {poly_deg}): R2_raw_poly={r2_raw_poly:.6g}\n"
            )

    # -------------------------------------------------------------------------
    # GNUPLOT scripts/plots (legacy)
    # -------------------------------------------------------------------------
    write_plot_target(results_dir, csv_path.name, f"target_exponents_{tag}Hz.png", target_freq,
                      source_label, n_vn, r2_RI_vn, p_vn, r2_ICI_vn)

    if vn_smin and vn_smin > 0:
        write_plot_target_smin(
            results_dir, csv_path.name, f"target_exponents_Smin_{tag}Hz.png", target_freq,
            source_label, vn_smin, n_vn_smin, r2_RI_vn_smin, p_vn_smin, r2_ICI_vn_smin
        )

    write_plot_index(results_dir, csv_path.name, f"index_exponents_{tag}Hz.png", target_freq,
                     source_label, n_idx, r2_ri, p_idx, r2_ici)

    write_plot_raw(results_dir, csv_path.name, f"exponents_raw_{tag}Hz.png", target_freq,
                   source_label, lin_rho, lin_im, log_fit_acc)

    write_plot_logistic_only(results_dir, csv_path.name, f"exponents_logistic_only_{tag}Hz.png", target_freq,
                             source_label, log_fit_acc)

    write_plot_logistic_diag(results_dir, csv_path.name, f"exponents_logistic_diagnostics_{tag}Hz.png",
                             target_freq, source_label, has_log=(log_fit_acc is not None))

    if poly_deg and poly_deg >= 2:
        write_plot_target_poly(
            results_dir, csv_path.name, f"target_exponents_polydeg{poly_deg}_{tag}Hz.png", target_freq,
            source_label, poly_deg, r2_RI_vn_poly, r2_ICI_vn_poly
        )
        if vn_smin and vn_smin > 0:
            write_plot_target_poly_smin(
                results_dir, csv_path.name, f"target_exponents_polydeg{poly_deg}_Smin_{tag}Hz.png", target_freq,
                source_label, poly_deg, vn_smin, r2_RI_vn_poly_smin, r2_ICI_vn_poly_smin
            )
        if poly_raw:
            write_plot_raw_poly(
                results_dir, csv_path.name, f"exponents_raw_polydeg{poly_deg}_{tag}Hz.png", target_freq,
                source_label, poly_deg, r2_raw_poly
            )

    # -------------------------------------------------------------------------
    # NEW: NONORM real-space plots (explicit names)
    # -------------------------------------------------------------------------
    # Decide x axis name for filenames
    x_tok = "THETA" if (x_col == "theta") else "S"

    # Sigma_imag plot
    out_sigma = f"NONORM_SIGMA_REAL_vs_{x_tok}_ALL_{tag}Hz.png"
    write_plot_nonorm_real(
        results_dir, csv_path.name, out_sigma, target_freq, source_label,
        x_col=x_col, y_col="sigma_imag",
        y_fit_cols={"LIN":"sigma_nonorm_fit_lin", "EXP":"sigma_nonorm_fit_exp", "POW":"sigma_nonorm_fit_pow"},
        r2_map={
            "LIN": (fit_nonorm_sigma["lin"]["r2"] if fit_nonorm_sigma.get("lin") else np.nan),
            "EXP": (fit_nonorm_sigma["exp"]["r2"] if fit_nonorm_sigma.get("exp") else np.nan),
            "POW": (fit_nonorm_sigma["pow"]["r2"] if fit_nonorm_sigma.get("pow") else np.nan),
        }
    )

    # Rho plot
    out_rho = f"NONORM_RHO_REAL_vs_{x_tok}_ALL_{tag}Hz.png"
    write_plot_nonorm_real(
        results_dir, csv_path.name, out_rho, target_freq, source_label,
        x_col=x_col, y_col="rho",
        y_fit_cols={"LIN":"rho_nonorm_fit_lin", "EXP":"rho_nonorm_fit_exp", "POW":"rho_nonorm_fit_pow"},
        r2_map={
            "LIN": (fit_nonorm_rho["lin"]["r2"] if fit_nonorm_rho.get("lin") else np.nan),
            "EXP": (fit_nonorm_rho["exp"]["r2"] if fit_nonorm_rho.get("exp") else np.nan),
            "POW": (fit_nonorm_rho["pow"]["r2"] if fit_nonorm_rho.get("pow") else np.nan),
        }
    )

    # Sigma_real plot (optional but written; may be mostly redundant)
    out_sigR = f"NONORM_SIGMAREAL_REAL_vs_{x_tok}_ALL_{tag}Hz.png"
    write_plot_nonorm_real(
        results_dir, csv_path.name, out_sigR, target_freq, source_label,
        x_col=x_col, y_col="sigma_real",
        y_fit_cols={"LIN":"sigma_real_nonorm_fit_lin", "EXP":"sigma_real_nonorm_fit_exp", "POW":"sigma_real_nonorm_fit_pow"},
        r2_map={
            "LIN": (fit_nonorm_sigR["lin"]["r2"] if fit_nonorm_sigR.get("lin") else np.nan),
            "EXP": (fit_nonorm_sigR["exp"]["r2"] if fit_nonorm_sigR.get("exp") else np.nan),
            "POW": (fit_nonorm_sigR["pow"]["r2"] if fit_nonorm_sigR.get("pow") else np.nan),
        }
    )

    # -------------------------------------------------------------------------
    # EXECUTE GNUPLOT (unless dry-run)
    # -------------------------------------------------------------------------
    if not dry_run:
        # legacy scripts
        _run_gnuplot("plot_target_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        if vn_smin and vn_smin > 0:
            _run_gnuplot("plot_target_exponents_smin.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_index_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_raw_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_logistic_only.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_logistic_diagnostics.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

        if poly_deg and poly_deg >= 2:
            _run_gnuplot("plot_target_exponents_poly.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
            if vn_smin and vn_smin > 0:
                _run_gnuplot("plot_target_exponents_poly_smin.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
            if poly_raw:
                _run_gnuplot("plot_raw_exponents_poly.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)

        # NEW NONORM scripts are named plot_<stem>.gp
        _run_gnuplot(f"plot_{Path(out_sigma).stem}.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)
        _run_gnuplot(f"plot_{Path(out_rho).stem}.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)
        _run_gnuplot(f"plot_{Path(out_sigR).stem}.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

    # -------------------------------------------------------------------------
    # NEW: PER-RUN EXCEL FULL AUDIT
    # -------------------------------------------------------------------------
    # Fit params table (one row per method-family)
    fit_rows = []

    def _add_fit_row(family, xname, yname, model, params_dict):
        row = {
            "family": family,
            "x": xname,
            "y": yname,
            "model": model,
        }
        if params_dict is None:
            row["ok"] = False
            fit_rows.append(row)
            return
        row["ok"] = True
        row.update(params_dict)
        fit_rows.append(row)

    # legacy summary rows (minimal, but traceable)
    _add_fit_row("RAW_LOG", "logS", "logRho", "linear", {"slope": slope_r, "intercept": intercept_r, "R2": r2_r, "p": p_r, "AIC": aic_lin_r})
    _add_fit_row("RAW_LOG", "logS", "logImag", "linear", {"slope": slope_i, "intercept": intercept_i, "R2": r2_i, "p": p_i, "AIC": aic_lin_i})
    _add_fit_row("INDEX_GOAL", "logS", "logRhoIdx", "through_origin", {"slope": a_ri, "n_idx": n_idx, "R2": r2_ri})
    _add_fit_row("INDEX_GOAL", "logS", "logImagIdx", "through_origin", {"slope": a_ici, "p_idx": p_idx, "R2": r2_ici})
    _add_fit_row("TARGET_VN", "logS", "logRI_vn", "linear", {"slope": slope_RI_vn, "intercept": int_RI_vn, "n_vn": n_vn, "R2": r2_RI_vn})
    _add_fit_row("TARGET_VN", "logS", "logICI_vn", "linear", {"slope": slope_ICI_vn, "intercept": int_ICI_vn, "p_vn": p_vn, "R2": r2_ICI_vn})

    if vn_smin and vn_smin > 0:
        _add_fit_row("TARGET_VN_SMIN", "logS", "logRI_vn", "linear", {"n_vn_smin": n_vn_smin, "R2": r2_RI_vn_smin})
        _add_fit_row("TARGET_VN_SMIN", "logS", "logICI_vn", "linear", {"p_vn_smin": p_vn_smin, "R2": r2_ICI_vn_smin})

    if poly_deg and poly_deg >= 2:
        _add_fit_row("TARGET_VN", "logS", "logRI_vn", f"poly_deg{poly_deg}", {"R2": r2_RI_vn_poly})
        _add_fit_row("TARGET_VN", "logS", "logICI_vn", f"poly_deg{poly_deg}", {"R2": r2_ICI_vn_poly})
        if vn_smin and vn_smin > 0:
            _add_fit_row("TARGET_VN_SMIN", "logS", "logRI_vn", f"poly_deg{poly_deg}", {"R2": r2_RI_vn_poly_smin})
            _add_fit_row("TARGET_VN_SMIN", "logS", "logICI_vn", f"poly_deg{poly_deg}", {"R2": r2_ICI_vn_poly_smin})
        if poly_raw:
            _add_fit_row("RAW_LOG", "logS", "logRho", f"poly_deg{poly_deg}", {"R2": r2_raw_poly})

    # logistic summary
    _add_fit_row("RAW_LOG", "logS", "logRho", "logistic4", (log_fit_acc if log_fit_acc is not None else None))

    # NEW: NONORM real-space fits (conductivity/moisture goal)
    def _add_nonorm_fitpack(yname, pack, prefix):
        xname = pack.get("x_col", x_col)
        _add_fit_row("NONORM_REAL", xname, yname, "linear", pack.get("lin"))
        _add_fit_row("NONORM_REAL", xname, yname, "exponential", pack.get("exp"))
        _add_fit_row("NONORM_REAL", xname, yname, "powerlaw", pack.get("pow"))

    _add_nonorm_fitpack("sigma_imag", fit_nonorm_sigma, "sigma_nonorm")
    _add_nonorm_fitpack("rho", fit_nonorm_rho, "rho_nonorm")
    _add_nonorm_fitpack("sigma_real", fit_nonorm_sigR, "sigma_real_nonorm")

    notes_rows = [{
        "run_root": str(run_root),
        "source": source_label,
        "sip_file": str(sip_xlsx),
        "logbook_file": str(log_xlsx),
        "logbook_sheet": used_sheet,
        "logbook_overlap": overlap_cnt,
        "logbook_block_rows": block_rows,
        "points_used": int(len(sip_df)),
        "x_for_nonorm": x_col,
        "ici_mode": ici_mode,
        "vn_smin": float(vn_smin) if (vn_smin and vn_smin > 0) else 0.0,
        "poly_deg": int(poly_deg) if (poly_deg and poly_deg >= 2) else 0,
        "poly_raw": bool(poly_raw),
        "logistic_status": log_status,
        "logistic_accepted": bool(log_fit_acc is not None),
        "logistic_reason": log_accept_reason,
    }]

    # Per-run Excel
    xlsx_path = results_dir / f"full_audit_{tag}Hz.xlsx"
    try:
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as w:
            sip_df.to_excel(w, sheet_name="data_long", index=False)
            pd.DataFrame(fit_rows).to_excel(w, sheet_name="fit_params", index=False)
            pd.DataFrame(notes_rows).to_excel(w, sheet_name="notes", index=False)
    except Exception as e:
        print(f"[WARN] failed to write per-run excel {xlsx_path}: {e}", flush=True)

    # Return summary + also return data for master workbook
    return {
        "run": str(run_root),
        "source": source_label,
        "status": "OK",
        "points": int(len(sip_df)),
        "logbook_sheet": used_sheet,
        "logbook_overlap": overlap_cnt,
        "logbook_block_rows": block_rows,
        "logistic_attempted": bool(log_fit is not None),
        "logistic_accepted": bool(log_fit_acc is not None),
        "logistic_status": log_status,
        "logistic_reason": log_accept_reason,
        "ici_mode": ici_mode,

        "n_vn": float(n_vn) if np.isfinite(n_vn) else n_vn,
        "p_vn": float(p_vn) if np.isfinite(p_vn) else p_vn,
        "n_idx": float(n_idx) if np.isfinite(n_idx) else n_idx,
        "p_idx": float(p_idx) if np.isfinite(p_idx) else p_idx,

        "x_nonorm": x_col,

        # nonorm key R2s (real space) for quick scanning
        "R2_nonorm_sigma_lin": (fit_nonorm_sigma["lin"]["r2"] if fit_nonorm_sigma.get("lin") else np.nan),
        "R2_nonorm_sigma_exp": (fit_nonorm_sigma["exp"]["r2"] if fit_nonorm_sigma.get("exp") else np.nan),
        "R2_nonorm_sigma_pow": (fit_nonorm_sigma["pow"]["r2"] if fit_nonorm_sigma.get("pow") else np.nan),

        "per_run_excel": str(xlsx_path),
        "joined_csv": str(csv_path),
        "_sip_df_for_master": sip_df,
        "_fit_rows_for_master": fit_rows,
    }


# =============================================================================
# MAIN CLI (UPDATED: write MASTER EXCEL for all runs)
# =============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("top", help="Top folder containing run subfolders")
    ap.add_argument("--freq", type=float, default=0.01, help="Target frequency in Hz")
    ap.add_argument("--sample", type=str, default=None,
                    help="Optional substring filter on SIP sheet names (for debugging)")
    ap.add_argument("--dry-run", action="store_true", help="Write scripts/CSVs but do not run gnuplot")
    ap.add_argument("--min-points", type=int, default=4, help="Minimum matched points to accept a run")
    ap.add_argument("--max-df", type=float, default=0.0, help="Max allowed |f_found - f_target| (0 disables)")
    ap.add_argument("--no-logistic", action="store_true", help="Disable logistic fit/plots")
    ap.add_argument("--logistic-min-points", type=int, default=6, help="Min points needed to attempt logistic")
    ap.add_argument("--logistic-fast", action="store_true", help="Use SciPy-free fast logistic search")
    ap.add_argument("--logistic-max-seconds", type=float, default=8.0, help="Time cap for logistic fit per run")
    ap.add_argument("--gnuplot-timeout", type=int, default=30, help="Timeout seconds per gnuplot call")
    ap.add_argument("--clean-results", action="store_true",
                    help="Delete files inside each run's Results/ before writing new outputs.")

    ap.add_argument("--ici-mode", type=str, default="sigref_over_sig",
                    choices=["sigref_over_sig","sig_over_sigref"],
                    help=("GOAL/INDEX ICI definition. TARGET is fixed to VN-style (sigref/sig) "
                          "so that both target panels trend downward."))

    ap.add_argument("--vn-smin", type=float, default=0.0,
                    help="Extra VN-style TARGET fit/plot using only points with S >= vn_smin.")

    ap.add_argument("--poly-deg", type=int, default=0,
                    help="Polynomial degree for extra fits/plots (>=2). 0 disables.")
    ap.add_argument("--poly-raw", action="store_true",
                    help="Also plot polynomial fit for RAW logRho vs logS (requires --poly-deg>=2).")

    args = ap.parse_args()

    base = Path(args.top).expanduser().resolve()
    if not base.exists():
        print(f"Not found: {base}")
        raise SystemExit(2)

    roots = discover_run_roots_fast(base)
    print(f"Found {len(roots)} run folder(s).", flush=True)

    rows = []
    all_data_long = []
    all_fit_params = []

    t0 = time.time()

    for i, run in enumerate(roots, 1):
        print(f"[{i}/{len(roots)}] START {run}", flush=True)
        ti = time.time()
        try:
            res = process_run(
                run, args.freq, args.sample, args.dry_run,
                args.min_points, args.max_df,
                do_logistic=(not args.no_logistic),
                logistic_min_points=args.logistic_min_points,
                logistic_fast=args.logistic_fast,
                logistic_max_seconds=args.logistic_max_seconds,
                gnuplot_timeout=args.gnuplot_timeout,
                clean_results=args.clean_results,
                ici_mode=args.ici_mode,
                vn_smin=args.vn_smin,
                poly_deg=args.poly_deg,
                poly_raw=args.poly_raw
            )
            if res is None:
                continue

            dt = time.time() - ti
            rows.append({k: v for k, v in res.items() if not k.startswith("_")})

            # Collect master tables
            sip_df = res.get("_sip_df_for_master", None)
            if isinstance(sip_df, pd.DataFrame):
                dfm = sip_df.copy()
                dfm.insert(0, "run", str(run))
                dfm.insert(1, "source", res.get("source", ""))
                all_data_long.append(dfm)

            fit_rows = res.get("_fit_rows_for_master", None)
            if isinstance(fit_rows, list) and fit_rows:
                for fr in fit_rows:
                    fr2 = dict(fr)
                    fr2["run"] = str(run)
                    fr2["source"] = res.get("source", "")
                    fr2["freq_hz"] = float(args.freq)
                    all_fit_params.append(fr2)

            msg = (f"[{i}/{len(roots)}] OK points={res.get('points',0)} "
                   f"n_vn={res.get('n_vn',np.nan):.3g} p_vn={res.get('p_vn',np.nan):.3g} "
                   f"x_nonorm={res.get('x_nonorm','')} "
                   f"R2sigma_lin={res.get('R2_nonorm_sigma_lin',np.nan):.3g} ({dt:.1f}s)")
            print(msg, flush=True)

        except subprocess.TimeoutExpired:
            rows.append({"run": str(run), "status": "FAIL", "reason": "gnuplot timeout", "points": 0})
            print(f"[{i}/{len(roots)}] FAIL gnuplot timeout", flush=True)
        except Exception as e:
            rows.append({"run": str(run), "status": "FAIL", "reason": str(e), "points": 0})
            print(f"[{i}/{len(roots)}] FAIL {e}", flush=True)

    # Summary CSV (legacy)
    tag = str(args.freq).replace(".","p")
    summary_path = base / f"Results_Summary_{tag}Hz.csv"
    if rows:
        pd.DataFrame(rows).to_csv(summary_path, index=False)
        print(f"\nSummary written: {summary_path}", flush=True)

    # NEW: MASTER EXCEL for all runs
    master_xlsx = base / f"Results_AllRuns_{tag}Hz.xlsx"
    try:
        with pd.ExcelWriter(master_xlsx, engine="openpyxl") as w:
            if all_data_long:
                pd.concat(all_data_long, ignore_index=True).to_excel(w, sheet_name="all_data_long", index=False)
            if all_fit_params:
                pd.DataFrame(all_fit_params).to_excel(w, sheet_name="all_fit_params", index=False)
            pd.DataFrame(rows).to_excel(w, sheet_name="summary", index=False)
        print(f"Master Excel written: {master_xlsx}", flush=True)
    except Exception as e:
        print(f"[WARN] failed to write master excel {master_xlsx}: {e}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)

if __name__ == "__main__":
    main()
