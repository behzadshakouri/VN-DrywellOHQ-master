#!/usr/bin/env python3
"""
auto_sip_pipeline_v27_noscipy.py  (FULLY DOCUMENTED / COMMENTED)

PURPOSE
-------
This script batch-processes SIP column experiment folders to compute:
  1) RAW exponent fits (logRho vs logS, logImag vs logS) + optional logistic (logRho)
  2) GOAL/INDEX through-origin fits on:
        logRI  = log10(rho/rho_ref)
        logICI = log10(ICI)
     where GOAL ICI can be defined as:
        - sigref_over_sig  (default): ICI = sigma_ref / sigma
        - sig_over_sigref            : ICI = sigma / sigma_ref
  3) TARGET "Van Nuys style" exponents in a way that matches the reference plots:
        - both panels trend downward
        - uses intercept (NOT forced through origin)
        - ALWAYS uses ICI = sigma_ref / sigma (downtrend by definition)
     i.e.
        logRI_vn  = logRho - logRho_ref
        logICI_vn = logSig_ref - logSig

In addition, v27 adds OPTIONAL extra fitting methods:
  - Polynomial regression (degree >=2) for TARGET VN curves (RI and ICI)
  - Polynomial regression for RAW logRho (optional)

The implementation avoids SciPy by including:
  - incomplete beta / continued fraction approximation (F-test p-value)
  - a fast grid-search logistic fitter

INPUT ASSUMPTIONS
-----------------
Directory structure per run:
  <run_root>/
     Analysis/         -> contains SIP output Excel (hint contains "output")
     SIP_Raw_Data/     -> contains logbook Excel (hint contains "Log Book")
     Results/          -> will be created and filled with csv + plots + report

SIP output Excel:
  - sheets correspond to "measurement IDs"
  - each sheet has columns for frequency, resistivity, imaginary conductivity
    (we attempt name matching; if fails, fallback indices are used)

Logbook Excel:
  - may have messy/incorrect column names
  - may contain multiple blocks/sheets; we select the *best contiguous block*
    based on overlap between logbook measurement IDs and SIP sheet names
  - saturation can come from:
      - a saturation column (direct)
      - OR computed from Mw, phi, Vcol (if saturation column absent)

MOST IMPORTANT FIX FOR "typos"
-----------------------------
We normalize measurement IDs so that "M 1", "M-1", "M_1", "m1" all match.
And we prefer the logbook block that overlaps the SIP sheet names the most.

OUTPUTS (per run folder)
------------------------
Results/joined_<tag>Hz.csv
Results/fit_report_<tag>Hz.txt
Results/target_exponents_<tag>Hz.png
Results/index_exponents_<tag>Hz.png
Results/exponents_raw_<tag>Hz.png
Results/exponents_logistic_only_<tag>Hz.png
Results/exponents_logistic_diagnostics_<tag>Hz.png

Optional outputs:
Results/target_exponents_Smin_<tag>Hz.png                       (if --vn-smin)
Results/target_exponents_polydegD_<tag>Hz.png                   (if --poly-deg D)
Results/target_exponents_polydegD_Smin_<tag>Hz.png              (if --poly-deg and --vn-smin)
Results/exponents_raw_polydegD_<tag>Hz.png                      (if --poly-raw with --poly-deg)

Also:
Results/debug_before_filter_<tag>Hz.csv
Results/dropped_rows_<tag>Hz.csv                                (if any invalid rows)

Summary output at base folder:
Results_Summary_<tag>Hz.csv

USAGE EXAMPLES
--------------
# baseline (same outputs as v27)
python3 auto_sip_pipeline_v27_noscipy.py "/path/to/SIP" \
  --freq 0.01 --min-points 4 \
  --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
  --clean-results --ici-mode sig_over_sigref

# add extra VN-style fit using S >= 0.30
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

# Expected folder names inside each run directory
ANALYSIS_DIR = "Analysis"
RAW_DIR      = "SIP_Raw_Data"
RESULTS_DIR  = "Results"

# Hints used to locate the two input Excel files
ANALYSIS_XLSX_HINT = "output"     # SIP result file often contains "output"
LOGBOOK_XLSX_HINT  = "Log Book"   # logbook file often contains "Log Book"

# Candidate logbook column names (we try to match, but also robust to errors)
LOGBOOK_ID_COLS_CAND     = ["Measurements","Measurement","Meas","MeasID","ID"]
LOGBOOK_SAT_COLS_CAND    = ["Degree Saturation M2","Degree Saturation","Saturation","Sw","S","DegreeSaturation"]
LOGBOOK_SAMPLE_COLS_CAND = ["Sample","Sample ID","SampleID","Soil Sample","Specimen"]

LOGBOOK_MASSWATER_CAND = ["Mass Water","Water Mass","Mass of Water","Mw","Water"]
LOGBOOK_POROSITY_CAND  = ["Soil Porosity","Porosity","phi","Φ"]
LOGBOOK_COLVOL_CAND    = ["Column Volume (m3)","Column Volume","Volume (m3)","Vcol"]

# Candidate SIP output names
SIP_FREQ_CANDS = ["Frequency","Frequency (Hz)","Freq","Hz"]
SIP_RHO_CANDS  = ["Resistivity","Resistivity (Ω·m)","Resistivity (Ohm m)","rho"]
SIP_IMAG_CANDS = ["Imaginary Conductivity","Imaginary Conductivity (µS/cm)","Imag","sigma_imag"]

# Logistic parameter bounds for stability in the grid-search
K_MIN = 0.05
K_MAX = 50.0

# Logistic acceptance logic (keep same as previous)
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
    Given a dataframe and a list of candidate column names, return the first
    matching column from the dataframe (case-insensitive). If none found, "".
    """
    low = {_norm(c): c for c in df.columns}
    for k in cands:
        kk = _norm(k)
        if kk in low:
            return low[kk]
    return ""

def find_first_matching_file(folder: Path, hint: str, ext=".xlsx"):
    """
    Find the first file under folder whose name contains 'hint' (case-insensitive)
    and ends with ext. Returns Path or None.
    """
    if not folder.exists():
        return None
    for p in folder.rglob(f"*{ext}"):
        if hint.lower() in p.name.lower():
            return p
    return None

def to_float(x):
    """
    Convert to float robustly:
      - handles strings including percent "40%"
      - returns NaN on failure
    """
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
    """
    Safe log10:
      - if x not finite or <=0 -> NaN
      - else log10(x)
    """
    x = float(x)
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return math.log10(x)

def sse(y, yhat):
    """Sum of squared errors."""
    r = y - yhat
    return float(np.sum(r*r))

def _is_blank(x) -> bool:
    """Return True if value is blank-ish (None, '', 'nan', 'none')."""
    if x is None:
        return True
    s = str(x).strip()
    return (s == "") or (s.lower() in ["nan", "none"])

def norm_measid(s: str) -> str:
    """
    Normalize measurement IDs (critical for matching SIP sheet names with logbook IDs).
    Example:
      "M 1" -> "m1"
      "M-1" -> "m1"
      "M_1" -> "m1"
    This collapses all separators so typos/mixed formatting don't break mapping.
    """
    s = str(s).strip().lower()
    s = re.sub(r'[^a-z0-9]+', '', s)
    return s

def pretty_source_name(run_folder_name: str) -> str:
    """
    Make a readable plot title from folder name:
      - underscores -> spaces
      - convert 'S-1' -> 'S_{-1}' for gnuplot enhanced text formatting
    """
    s = str(run_folder_name).replace("_", " ").strip()
    s = re.sub(r'\bS-(\d+)\b', r'S_{-\1}', s)
    return s


# =============================================================================
# F-DISTRIBUTION + P-VALUE (NO SCIPY)
# =============================================================================
# Used only to compute model p-values for linear/logistic fits.
# Implementation: incomplete beta function via continued fraction.

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
    Compute p-value for regression model using F-test.
    Here df_model is number of predictors (NOT counting intercept).
    For linear y=a*x+b => df_model=1.
    For logistic 4-params we approximate df_model=3 (effective).
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
    """
    Gaussian AIC:
      AIC = n*ln(ss_res/n) + 2*k_params
    """
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)


# =============================================================================
# FIT ROUTINES: LINEAR / THROUGH ORIGIN / POLYNOMIAL
# =============================================================================

def linear_fit(x, y):
    """
    Ordinary least squares y = slope*x + intercept.
    Returns:
      slope, intercept, r2, pvalue, ss_res, aic
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
    """
    Least squares forced through origin: y = a*x.
    Used for GOAL plots (your "goal" through-origin convention).
    Returns:
      a, r2, ss_res
    """
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
# LOGISTIC FIT (NO SCIPY)
# =============================================================================

def logistic4(x, L, k, x0, c):
    """
    4-parameter logistic:
       f(x)=c + L/(1+exp(-k*(x-x0)))
    """
    z = k*(x - x0)
    z = np.clip(z, -60.0, 60.0)  # avoid overflow in exp
    return c + L/(1.0 + np.exp(-z))

def _logistic_init_guess(x, y):
    """
    Heuristic initialization for logistic:
      - c ~ min(y)
      - L ~ max(y)-min(y)
      - x0 ~ x at y ~ mid
      - k from slope near mid (clipped)
    """
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
    "Fast" logistic fit using a coarse coordinate grid search around an initial guess.
    - This is not as accurate as a proper non-linear optimizer, but is robust and SciPy-free.
    - We stop early if time exceeds max_seconds.

    Returns:
      (fit_dict, status_str)
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

    # step sizes around guess; shrink each iteration
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

    # coordinate search loop
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

        # shrink steps
        stepL *= 0.6; stepk *= 0.6; stepx0 *= 0.6; stepc *= 0.6
        if not improved and max(stepL, stepk, stepx0, stepc) < 1e-4:
            break

    # compute metrics
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
    """
    Decide whether to accept logistic fit over linear.
    Criteria:
      - if R2 improves by LOGISTIC_R2_EPS
      - OR (optional) if AIC improves
    """
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
      - loops over all sheets (each sheet = measurement)
      - finds the row closest to target_freq
      - extracts rho and sigma_imag at that frequency
      - returns dataframe with:
          measurement, measurement_norm, freq_hz, rho, sigma_imag

    sheet_filter:
      if provided, only sheets containing this substring are processed.
      (useful for debugging a single run).
    """
    xls = pd.ExcelFile(sip_xlsx)
    out = []
    for sh in xls.sheet_names:
        if sheet_filter and (sheet_filter.lower() not in str(sh).lower()):
            continue

        df = pd.read_excel(sip_xlsx, sheet_name=sh)
        if df.empty:
            continue

        # Attempt to identify key columns by name; fallback if needed
        freq_col = pick_col(df, SIP_FREQ_CANDS) or df.columns[0]
        rho_col  = pick_col(df, SIP_RHO_CANDS)
        imag_col = pick_col(df, SIP_IMAG_CANDS)

        # If named columns not found, fallback to legacy indices (as in earlier versions)
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

        # If user provides max_df, skip sheets too far away from target frequency
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
# LOGBOOK READING: ROBUST BLOCK SELECTION
# =============================================================================
# Key idea:
#   Logbooks often contain repeated measurement blocks, notes, empty rows, etc.
#   We split each sheet into contiguous blocks of rows where both:
#     - measurement ID exists
#     - saturation info exists (either direct S column OR Mw/phi/V present)
#   Then we choose the block maximizing overlap with SIP measurement IDs.

def _contiguous_blocks(mask: np.ndarray):
    """
    Given a boolean array, return list of contiguous (i0,i1) index ranges where mask is True.
    """
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
    Read logbook and choose best sheet/block by overlap with SIP measurement IDs.

    Returns:
      sat_map: dict measurement_norm -> S (float)
      used_sheet: sheet name
      overlap_cnt: number of IDs overlapping SIP sheets
      block_rows: number of rows in selected block
    """
    xls = pd.ExcelFile(log_xlsx)
    sip_set = set([m for m in sip_measurements_norm if not _is_blank(m)])

    best = None
    # best tuple: (score, sheet, blk_df, id_col, sat_col, mw_col, por_col, vol_col, overlap_count)

    for sh in xls.sheet_names:
        try:
            df = pd.read_excel(log_xlsx, sheet_name=sh)
        except Exception:
            continue
        if df is None or df.empty:
            continue

        # Identify measurement ID column (required)
        id_col  = pick_col(df, LOGBOOK_ID_COLS_CAND) or ""

        # Saturation may be direct column OR computed from Mw, phi, Vcol
        sat_col = pick_col(df, LOGBOOK_SAT_COLS_CAND)
        mw_col  = pick_col(df, LOGBOOK_MASSWATER_CAND)
        por_col = pick_col(df, LOGBOOK_POROSITY_CAND)
        vol_col = pick_col(df, LOGBOOK_COLVOL_CAND)

        if not id_col:
            continue

        has_meas = ~df[id_col].astype(str).apply(_is_blank)

        # Determine if row has saturation information
        if sat_col:
            has_sat = df[sat_col].apply(lambda x: np.isfinite(to_float(x)))
        else:
            # require Mw+phi+V to compute S
            if not (mw_col and por_col and vol_col):
                continue
            has_sat = df.apply(
                lambda r: np.isfinite(to_float(r.get(mw_col, np.nan)))
                          and np.isfinite(to_float(r.get(por_col, np.nan)))
                          and np.isfinite(to_float(r.get(vol_col, np.nan))),
                axis=1
            )

        usable = (has_meas & has_sat).to_numpy(bool)

        # Split into contiguous usable blocks
        blocks = _contiguous_blocks(usable)
        if not blocks:
            continue

        # Score each block by overlap with SIP IDs
        for (i0, i1) in blocks:
            blk = df.iloc[i0:i1+1].copy()
            blk_ids = blk[id_col].astype(str).apply(norm_measid).tolist()
            overlap = len(set(blk_ids) & sip_set)

            # scoring: prioritize overlap heavily
            score = overlap * 1000 + (i1 - i0 + 1)
            if (best is None) or (score > best[0]):
                best = (score, sh, blk, id_col, sat_col, mw_col, por_col, vol_col, overlap)

    if best is None:
        raise RuntimeError(f"Logbook has no usable saturation blocks. Sheets: {xls.sheet_names}")

    _, sh, blk, id_col, sat_col, mw_col, por_col, vol_col, overlap = best

    # Build saturation map (normalized measurement ID -> S)
    sat_map = {}
    for _, r in blk.iterrows():
        mid = str(r.get(id_col, "")).strip()
        if _is_blank(mid):
            continue

        S = np.nan
        if sat_col:
            # S might be in fraction or percent
            S = to_float(r.get(sat_col, np.nan))
            if np.isfinite(S) and S > 1.0 and S <= 100.0:
                S = S/100.0
        else:
            # compute S = Mw / (rho_w * phi * V)
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
    """
    Delete files in Results/ (do not delete folder).
    This ensures old plots/csv do not remain.
    """
    if not results_dir.exists():
        return
    for p in results_dir.iterdir():
        try:
            if p.is_file() or p.is_symlink():
                p.unlink()
        except Exception:
            pass


# =============================================================================
# GNUPLOT SCRIPT WRITERS
# =============================================================================
# These produce .gp scripts (using column("name") so headers are robust),
# then gnuplot is executed to create PNG images.

def write_plot_target(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                      source_label: str,
                      n_vn: float, r2_ri_vn: float, p_vn: float, r2_ici_vn: float):
    """
    TARGET plot matching Van Nuys style:
      - Two panels
      - both downward trend
      - uses VN-defined curves with intercept
      - titles use folder-based source label
    """
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
    """
    Optional: TARGET plot but only points with S>=vn_smin.
    Helpful when low-S points are wrong/noisy.
    """
    label = "{:g}".format(target_freq)

    # If not enough points, write a placeholder plot to avoid gnuplot errors.
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
    """
    Optional: polynomial fit overlay on TARGET curves (full set).
    """
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
    """
    Optional: polynomial fit overlay but only for S>=vn_smin.
    """
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
    """
    GOAL plot: through-origin fit on logRhoIdx/logImagIdx vs logS.
    This uses --ici-mode for ICI definition.
    """
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
    """
    RAW plot:
      - panel 1: logRho vs logS with linear + optional logistic overlay
      - panel 2: logImag vs logS with linear
    """
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
    """
    Optional: polynomial fit overlay on RAW logRho vs logS.
    """
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
    """
    Logistic-only plot:
      - always plots data
      - if logistic accepted, overlay logistic curve
      - else the curve is suppressed
    """
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
    """
    Logistic diagnostics plot (3 panels):
      1) logRho vs logS with linear + logistic (if accepted)
      2) residuals of linear
      3) residuals of logistic (or a placeholder if not accepted)
    """
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
# DISCOVERY + GNUPLOT EXECUTION
# =============================================================================

def discover_run_roots_fast(base: Path):
    """
    Find run roots by locating folders that contain:
      - Analysis/ directory
      - and a sibling SIP_Raw_Data/ directory
    """
    roots = set()
    for a in base.rglob(ANALYSIS_DIR):
        if not a.is_dir():
            continue
        run = a.parent
        if (run / RAW_DIR).is_dir():
            roots.add(run)
    return sorted(roots)

def _run_gnuplot(script_name: str, cwd: Path, timeout_s: int, fatal: bool):
    """
    Execute gnuplot script in cwd.
    If fatal=True, raise error; else print warning and continue.
    """
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
# MAIN PROCESSING FOR ONE RUN
# =============================================================================

def process_run(run_root: Path, target_freq: float, sheet_filter: str|None,
                dry_run: bool, min_points: int, max_df: float,
                do_logistic: bool, logistic_min_points: int,
                logistic_fast: bool, logistic_max_seconds: float,
                gnuplot_timeout: int,
                clean_results: bool,
                ici_mode: str,
                vn_smin: float,
                poly_deg: int,
                poly_raw: bool):
    """
    Process a single run folder:
      - locate SIP output + logbook
      - read SIP data at target freq
      - choose best logbook block and map S to SIP measurement IDs
      - compute logs and fits
      - write outputs: joined CSV, report, plots

    Returns a dictionary row for summary CSV.
    """
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    # Locate required Excel files
    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return None

    # Source label is run folder name (requested by you)
    source_label = pretty_source_name(run_root.name)

    # Ensure Results/ exists and optionally clean
    results_dir.mkdir(parents=True, exist_ok=True)
    if clean_results:
        clean_results_dir(results_dir)

    # Read SIP sheets
    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df, sheet_filter=sheet_filter)
    tag = str(target_freq).replace(".","p")

    if sip_df.empty:
        return {"run": str(run_root), "status": "FAIL",
                "reason": "No SIP sheets readable / target freq missing",
                "points": 0}

    # Read logbook; select best block by overlap with SIP measurement IDs
    sat_map, used_sheet, overlap_cnt, block_rows = read_logbook_bestblock_by_overlap(
        log_xlsx, sip_df["measurement_norm"].tolist()
    )

    # Map saturations onto SIP measurement rows
    sip_df["S"] = sip_df["measurement_norm"].map(sat_map).astype(float)

    # Compute logs
    sip_df["logS"]    = sip_df["S"].apply(safe_log10)
    sip_df["logRho"]  = sip_df["rho"].apply(safe_log10)
    sip_df["logImag"] = sip_df["sigma_imag"].apply(safe_log10)

    # Debug dump before filtering (helps find bad S/rho/sigma)
    dbg_all = sip_df.copy()
    dbg_all.to_csv(results_dir / f"debug_before_filter_{tag}Hz.csv", index=False)

    # Filter out invalid rows (log inputs must be finite)
    mask = np.isfinite(sip_df["logS"]) & np.isfinite(sip_df["logRho"]) & np.isfinite(sip_df["logImag"])
    dropped = sip_df.loc[~mask, ["measurement", "S", "rho", "sigma_imag"]].copy()
    if len(dropped) > 0:
        dropped.to_csv(results_dir / f"dropped_rows_{tag}Hz.csv", index=False)
    sip_df = sip_df.loc[mask].copy()

    # Enforce minimum points after filtering
    if len(sip_df) < min_points:
        return {"run": str(run_root), "status": "FAIL",
                "reason": f"Not enough matched points (need >= {min_points})",
                "points": int(len(sip_df)),
                "logbook_sheet": used_sheet,
                "logbook_overlap": overlap_cnt,
                "logbook_block_rows": block_rows,
                "source": source_label}

    # Choose reference point as max S (your convention)
    iref = int(np.nanargmax(sip_df["S"].to_numpy(float)))
    rho_ref  = float(sip_df.iloc[iref]["rho"])
    sig_ref  = float(sip_df.iloc[iref]["sigma_imag"])
    logRho_ref = safe_log10(rho_ref)
    logSig_ref = safe_log10(sig_ref)

    # GOAL/INDEX definitions
    sip_df["rho_idx"]   = sip_df["rho"] / rho_ref
    sip_df["logRhoIdx"] = sip_df["logRho"] - logRho_ref

    if ici_mode == "sig_over_sigref":
        # ICI = sig/sig_ref => logICI = logSig - logSig_ref
        sip_df["imag_idx"]   = sip_df["sigma_imag"] / sig_ref
        sip_df["logImagIdx"] = sip_df["logImag"] - logSig_ref
    else:
        # ICI = sig_ref/sig => logICI = logSig_ref - logSig
        sip_df["imag_idx"]   = sig_ref / sip_df["sigma_imag"]
        sip_df["logImagIdx"] = logSig_ref - sip_df["logImag"]

    # RAW fits (linear)
    slope_r, intercept_r, r2_r, p_r, ss_res_lin_r, aic_lin_r = linear_fit(sip_df["logS"], sip_df["logRho"])
    slope_i, intercept_i, r2_i, p_i, ss_res_lin_i, aic_lin_i = linear_fit(sip_df["logS"], sip_df["logImag"])
    lin_rho = {"slope": slope_r, "intercept": intercept_r, "r2": r2_r, "p": p_r, "aic": aic_lin_r}
    lin_im  = {"slope": slope_i, "intercept": intercept_i, "r2": r2_i, "p": p_i, "aic": aic_lin_i}

    # GOAL/INDEX through-origin fits
    a_ri,  r2_ri,  _ = linear_fit_through_origin(sip_df["logS"], sip_df["logRhoIdx"])
    a_ici, r2_ici, _ = linear_fit_through_origin(sip_df["logS"], sip_df["logImagIdx"])

    # Convert slopes to positive exponents (your convention)
    n_idx = -a_ri
    if ici_mode == "sig_over_sigref":
        # if logImagIdx = logSig - logSig_ref, slope is positive if sigma increases with S
        p_idx = +a_ici
    else:
        p_idx = -a_ici

    # TARGET (Van Nuys style): ALWAYS downtrend
    # This is independent of --ici-mode by design (so target always matches VN template).
    sip_df["logRI_vn"]  = sip_df["logRhoIdx"]
    sip_df["logICI_vn"] = logSig_ref - sip_df["logImag"]

    slope_RI_vn,  int_RI_vn,  r2_RI_vn,  _, _, _ = linear_fit(sip_df["logS"], sip_df["logRI_vn"])
    slope_ICI_vn, int_ICI_vn, r2_ICI_vn, _, _, _ = linear_fit(sip_df["logS"], sip_df["logICI_vn"])

    n_vn = -slope_RI_vn
    p_vn = -slope_ICI_vn

    sip_df["logRI_vn_fit"]  = slope_RI_vn  * sip_df["logS"] + int_RI_vn
    sip_df["logICI_vn_fit"] = slope_ICI_vn * sip_df["logS"] + int_ICI_vn

    # Optional TARGET fit using subset S>=vn_smin
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

    # Polynomial regression on TARGET VN curves
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

    # Polynomial regression on RAW logRho (optional)
    r2_raw_poly = np.nan
    sip_df["logRho_polyfit"] = np.nan
    if poly_raw and (poly_deg and poly_deg >= 2):
        cR, _, r2p, _ = poly_fit_predict(sip_df["logS"], sip_df["logRho"], poly_deg)
        r2_raw_poly = r2p
        if cR is not None:
            sip_df["logRho_polyfit"] = np.polyval(cR, sip_df["logS"].to_numpy(float))

    # Logistic fit on RAW logRho (optional)
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

    # Residuals and fitted curves for diagnostics plots
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

    # Write joined CSV with all important columns
    csv_path = results_dir / f"joined_{tag}Hz.csv"
    sip_df[[
        "measurement","measurement_norm","freq_hz","S","rho","sigma_imag",
        "logS","logRho","logImag",
        "rho_idx","imag_idx","logRhoIdx","logImagIdx",

        # TARGET VN columns
        "logRI_vn","logICI_vn","logRI_vn_fit","logICI_vn_fit",
        "logRI_vn_fit_smin","logICI_vn_fit_smin",

        # TARGET polynomial
        "logRI_vn_polyfit","logICI_vn_polyfit",
        "logRI_vn_polyfit_smin","logICI_vn_polyfit_smin",

        # RAW polynomial
        "logRho_polyfit",

        # diagnostics
        "resid_lin_rho","resid_log_rho",
        "logRho_linfit","logRho_logfit",
    ]].to_csv(csv_path, index=False)

    # Write report for traceability
    report = results_dir / f"fit_report_{tag}Hz.txt"
    report.write_text(
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
        f"Logistic: status={log_status}, accepted={bool(log_fit_acc is not None)}, reason={log_accept_reason}\n"
    )

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

    # Write gnuplot scripts (base)
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

    # Write gnuplot scripts (optional polynomial)
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

    # Execute gnuplot unless dry-run
    if not dry_run:
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

    # Return summary row
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

        "vn_smin": float(vn_smin) if (vn_smin and vn_smin > 0) else 0.0,
        "n_vn_smin": float(n_vn_smin) if np.isfinite(n_vn_smin) else n_vn_smin,
        "p_vn_smin": float(p_vn_smin) if np.isfinite(p_vn_smin) else p_vn_smin,

        "poly_deg": int(poly_deg) if (poly_deg and poly_deg >= 2) else 0,
        "R2_RI_vn_poly": float(r2_RI_vn_poly) if np.isfinite(r2_RI_vn_poly) else r2_RI_vn_poly,
        "R2_ICI_vn_poly": float(r2_ICI_vn_poly) if np.isfinite(r2_ICI_vn_poly) else r2_ICI_vn_poly,
        "R2_RI_vn_poly_smin": float(r2_RI_vn_poly_smin) if np.isfinite(r2_RI_vn_poly_smin) else r2_RI_vn_poly_smin,
        "R2_ICI_vn_poly_smin": float(r2_ICI_vn_poly_smin) if np.isfinite(r2_ICI_vn_poly_smin) else r2_ICI_vn_poly_smin,
        "poly_raw": bool(poly_raw),
        "R2_raw_poly": float(r2_raw_poly) if np.isfinite(r2_raw_poly) else r2_raw_poly,
    }


# =============================================================================
# MAIN CLI
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

    # Find all run folders
    roots = discover_run_roots_fast(base)
    print(f"Found {len(roots)} run folder(s).", flush=True)

    rows = []
    t0 = time.time()

    # Process each run folder
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
            rows.append(res)
            dt = time.time() - ti

            # Print compact per-run summary
            msg = (f"[{i}/{len(roots)}] OK points={res.get('points',0)} "
                   f"n_vn={res.get('n_vn',np.nan):.3g} p_vn={res.get('p_vn',np.nan):.3g}")

            if args.vn_smin and args.vn_smin > 0:
                msg += (f" | n_vn_smin={res.get('n_vn_smin',np.nan):.3g} "
                        f"p_vn_smin={res.get('p_vn_smin',np.nan):.3g}")

            if args.poly_deg and args.poly_deg >= 2:
                msg += (f" | poly{args.poly_deg} R2(RI)={res.get('R2_RI_vn_poly',np.nan):.3g} "
                        f"R2(ICI)={res.get('R2_ICI_vn_poly',np.nan):.3g}")
                if args.vn_smin and args.vn_smin > 0:
                    msg += (f" | poly{args.poly_deg}@Smin R2(RI)={res.get('R2_RI_vn_poly_smin',np.nan):.3g} "
                            f"R2(ICI)={res.get('R2_ICI_vn_poly_smin',np.nan):.3g}")
                if args.poly_raw:
                    msg += f" | rawpoly R2={res.get('R2_raw_poly',np.nan):.3g}"

            msg += (f" logistic_attempted={res.get('logistic_attempted',False)} "
                    f"logistic_accepted={res.get('logistic_accepted',False)} ({dt:.1f}s)")
            print(msg, flush=True)

        except subprocess.TimeoutExpired:
            rows.append({"run": str(run), "status": "FAIL", "reason": "gnuplot timeout", "points": 0})
            print(f"[{i}/{len(roots)}] FAIL gnuplot timeout", flush=True)
        except Exception as e:
            rows.append({"run": str(run), "status": "FAIL", "reason": str(e), "points": 0})
            print(f"[{i}/{len(roots)}] FAIL {e}", flush=True)

    # Write overall summary CSV at the base directory
    tag = str(args.freq).replace(".","p")
    summary_path = base / f"Results_Summary_{tag}Hz.csv"
    if rows:
        pd.DataFrame(rows).to_csv(summary_path, index=False)
        print(f"\nSummary written: {summary_path}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)

if __name__ == "__main__":
    main()
