#!/usr/bin/env python3
"""
auto_sip_pipeline_v27_noscipy.py  (UPDATED: normalize switch + real-space fits/plots + Excel dumps)

WHAT'S NEW IN THIS UPDATE
-------------------------
You asked to KEEP all previous approaches/plots and ADD:

1) A global switch to run EVERYTHING in:
   - NORMALIZED mode (default): uses a reference measurement (max S) to form indices like rho/rho_ref
   - UNNORMALIZED mode: does NOT divide by a reference measurement for "index/target" style quantities.

   CLI:
     --normalize          (default ON)
     --no-normalize       (turn OFF normalization)

   Notes on definitions when --no-normalize is used:
     - "RAW" fits are always inherently unnormalized (they use rho and sigma_imag directly)
       but are labeled consistently as "unorm" when normalization is off.
     - "INDEX/GOAL" (through-origin) and "TARGET/VN-style" require ref-based definitions
       in the classic form. When normalization is OFF, we compute and plot ALTERNATIVE
       "un-normalized" analogs that do not use rho_ref or sigma_ref:
         * RI_unorm  = rho                 (logRI_unorm = logRho)
         * ICI_unorm = 1/sigma_imag        (logICI_unorm = -logSigma)
       This keeps an "ICI-like" downtrend with increasing saturation without using a reference.

2) Write EVERYTHING to Excel:
   - Per folder/run:
       Results/all_values_<tag>Hz_<normtag>.xlsx
     Contains:
       * Data (real + log values, indices, VN/alt columns, predictions, residuals, etc.)
       * FitSummary (all fit parameters and metrics)
       * Meta (file names, chosen logbook block, overlap, etc.)
   - Global:
       <top>/Results_All_<tag>Hz_<normtag>.xlsx
     Contains:
       * Summary table across runs
       * (Optional) concatenated fit summaries per run

3) Add REAL-SPACE (no log) fits and plots:
   - linear fits on real values:
       rho vs S, sigma_imag vs S
       sigma_imag vs theta   (if theta can be computed)
   - optional polynomial fits (deg>=2) on real values
   - optional logistic fits on real-space rho vs S  (same SciPy-free fitter, just with x=S)
     (logistic on sigma is disabled by default; you can enable easily if you want)

4) Add fits WITHOUT normalization (as above) and produce lots of plots,
   with names that clearly encode:
     - norm/unorm
     - domain: loglog or real
     - which response: rho, sigma, RI, ICI, etc.
     - method: linear, polydegD, logistic

FINAL SCIENTIFIC GOAL
---------------------
You said the final goal is to find relation between electric conductivity and moisture content.
This script now explicitly produces:
  - sigma_imag vs S  (real and log-log)
  - if porosity phi is available per measurement: theta = S * phi, then sigma_imag vs theta
and exports those columns and fit parameters to Excel for each run and globally.

INPUT ASSUMPTIONS (same as v27)
-------------------------------
Per run:
  <run_root>/
     Analysis/         -> SIP output Excel (name contains "output")
     SIP_Raw_Data/     -> logbook Excel (name contains "Log Book")
     Results/          -> outputs

SIP output Excel sheets:
  - each sheet is a measurement ID

Logbook:
  - may have typos/messy columns
  - we choose the best contiguous block by overlap with SIP sheet names
  - saturation comes from either a saturation column or computed from Mw, phi, Vcol

USAGE EXAMPLES
--------------
# NORMALIZED (default), like your current runner:
python3 auto_sip_pipeline_v27_noscipy.py "/path/to/SIP" \
  --freq 0.01 --min-points 4 \
  --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
  --clean-results --ici-mode sig_over_sigref \
  --poly-deg 2 --poly-raw

# UNNORMALIZED mode (no ref division on goal/target; use RI=rho and ICI=1/sigma):
python3 auto_sip_pipeline_v27_noscipy.py "/path/to/SIP" \
  --freq 0.01 --min-points 4 \
  --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
  --clean-results --ici-mode sig_over_sigref \
  --poly-deg 2 --poly-raw \
  --no-normalize

# Also do real-space polynomial fits & plots (deg 2) and real-space logistic rho vs S:
python3 auto_sip_pipeline_v27_noscipy.py "/path/to/SIP" \
  --freq 0.01 --min-points 4 \
  --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
  --clean-results --ici-mode sig_over_sigref \
  --poly-deg 2 --poly-raw \
  --real-fits --real-logistic
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
SIP_IMAG_CANDS = ["Imaginary Conductivity","Imaginary Conductivity (µS/cm)","Imag","sigma_imag","sigma","conductivity"]

# Logistic bounds for stability
K_MIN = 0.05
K_MAX = 50.0

LOGISTIC_R2_EPS = 0.005
LOGISTIC_REQUIRE_BETTER_AIC = True

# =============================================================================
# BASIC HELPERS
# =============================================================================

def _norm(s: str) -> str:
    return str(s).strip().lower()

def pick_col(df: pd.DataFrame, cands):
    low = {_norm(c): c for c in df.columns}
    for k in cands:
        kk = _norm(k)
        if kk in low:
            return low[kk]
    return ""

def find_first_matching_file(folder: Path, hint: str, ext=".xlsx"):
    if not folder.exists():
        return None
    for p in folder.rglob(f"*{ext}"):
        if hint.lower() in p.name.lower():
            return p
    return None

def to_float(x):
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
    x = float(x)
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return math.log10(x)

def sse(y, yhat):
    r = y - yhat
    return float(np.sum(r*r))

def _is_blank(x) -> bool:
    if x is None:
        return True
    s = str(x).strip()
    return (s == "") or (s.lower() in ["nan", "none"])

def norm_measid(s: str) -> str:
    s = str(s).strip().lower()
    s = re.sub(r'[^a-z0-9]+', '', s)
    return s

def pretty_source_name(run_folder_name: str) -> str:
    s = str(run_folder_name).replace("_", " ").strip()
    # gnuplot enhanced-text: S-1 -> S_{-1}
    s = re.sub(r'\bS-(\d+)\b', r'S_{-\1}', s)
    return s

# =============================================================================
# F-DISTRIBUTION + P-VALUE (NO SCIPY)
# =============================================================================

def betacf(a, b, x, maxit=200, eps=3e-14):
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
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    ln_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a+b)
    bt = math.exp(math.log(x)*a + math.log(1.0-x)*b - ln_beta)
    if x < (a+1.0)/(a+b+2.0):
        return bt*betacf(a,b,x)/a
    else:
        return 1.0 - bt*betacf(b,a,1.0-x)/b

def f_dist_sf(F, d1, d2):
    if not np.isfinite(F) or F < 0:
        return np.nan
    x = (d1*F)/(d1*F + d2)
    cdf = betai(d1/2.0, d2/2.0, x)
    return max(0.0, 1.0 - cdf)

def f_test_pvalue(ss_res, ss_tot, df_model, n):
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
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)

# =============================================================================
# FIT ROUTINES
# =============================================================================

def linear_fit(x, y):
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
# LOGISTIC FIT (NO SCIPY) - works for x=logS or x=S equally
# =============================================================================

def logistic4(x, L, k, x0, c):
    z = k*(x - x0)
    z = np.clip(z, -60.0, 60.0)
    return c + L/(1.0 + np.exp(-z))

def _logistic_init_guess(x, y):
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

        # legacy fallback
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
# LOGBOOK READING: ROBUST BLOCK SELECTION + ALSO RETURN phi IF AVAILABLE
# =============================================================================

def _contiguous_blocks(mask: np.ndarray):
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
    Returns:
      maps: dict with keys:
        'S'   : measurement_norm -> S
        'phi' : measurement_norm -> phi  (if available)
      meta: used_sheet, overlap_cnt, block_rows
    """
    xls = pd.ExcelFile(log_xlsx)
    sip_set = set([m for m in sip_measurements_norm if not _is_blank(m)])

    best = None
    # best tuple: (score, sheet, blk_df, id_col, sat_col, mw_col, por_col, vol_col, overlap)

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

    S_map = {}
    phi_map = {}

    for _, r in blk.iterrows():
        mid = str(r.get(id_col, "")).strip()
        if _is_blank(mid):
            continue
        midn = norm_measid(mid)

        # phi (porosity) if present
        if por_col:
            ph = to_float(r.get(por_col, np.nan))
            if np.isfinite(ph) and ph > 0 and ph <= 1.5:
                phi_map[midn] = float(ph)

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
            S_map[midn] = float(S)

    if not S_map:
        raise RuntimeError(f"Selected logbook sheet '{sh}' but mapped 0 saturation rows.")

    return {"S": S_map, "phi": phi_map}, sh, int(overlap), int(len(blk))

# =============================================================================
# RESULTS CLEANING
# =============================================================================

def clean_results_dir(results_dir: Path):
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
# GNUPLOT SCRIPT WRITERS (LOG-LOG; KEEP YOUR EXISTING LOOK)
# =============================================================================
# We keep your v27 scripts as-is style-wise, but add norm tags and new plot types.

def _gp_term():
    return "set term pngcairo size 1200,1600 enhanced font 'Arial,26'\n"

def write_plot_target_loglog(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                             source_label: str, n_val: float, r2a: float, p_val: float, r2b: float,
                             col_y1: str, col_y1_fit: str, col_y2: str, col_y2_fit: str,
                             title1: str, title2: str, ylab1: str, ylab2: str,
                             legend1: str, legend2: str):
    """
    Generic 2-panel plot for "target-like" quantities in log-log domain:
      x = logS, y = log10(quantity)
    """
    label = "{:g}".format(target_freq)
    gp = f"""reset
set datafile separator ','
set datafile missing ''
{_gp_term()}set output '{out_png}'
set multiplot layout 2,1

set grid
set tics out
set key top right
set xlabel 'log10(Saturation)'

set title '{source_label}: {title1} at {label} Hz'
set ylabel '{ylab1}'
plot '{csv_name}' using (column("logS")):(column("{col_y1}")) with points pt 7 ps 1.8 title 'data', \\
     '{csv_name}' using (column("logS")):(column("{col_y1_fit}")) with lines dt 2 lw 3 title sprintf('{legend1}', {n_val}, {r2a})

set title '{source_label}: {title2} at {label} Hz'
set ylabel '{ylab2}'
plot '{csv_name}' using (column("logS")):(column("{col_y2}")) with points pt 7 ps 1.8 title 'data', \\
     '{csv_name}' using (column("logS")):(column("{col_y2_fit}")) with lines dt 2 lw 3 title sprintf('{legend2}', {p_val}, {r2b})

unset multiplot
"""
    (results_dir / "plot_target_loglog.gp").write_text(gp)

def write_plot_raw_loglog(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                          source_label: str,
                          lin_rho, lin_sig, log_fit_acc,
                          normtag: str):
    """
    RAW log-log plot:
      panel1: logRho vs logS with linear + (optional) logistic
      panel2: logSigma vs logS with linear
    """
    label = "{:g}".format(target_freq)
    a_r = lin_rho["slope"]; b_r = lin_rho["intercept"]; r2_r = lin_rho["r2"]; p_r = lin_rho["p"]
    a_s = lin_sig["slope"]; b_s = lin_sig["intercept"]; r2_s = lin_sig["r2"]; p_s = lin_sig["p"]

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

    gp = f"""reset
set datafile separator ','
set datafile missing ''
{_gp_term()}set output '{out_png}'
set multiplot layout 2,1

set grid
set key top right

set title '{source_label}: RAW (log-log, {normtag}) - log10(rho) vs log10(S) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
a1={a_r:.12g}
b1={b_r:.12g}
f1(x)=a1*x+b1
{log_defs}
plot '{csv_name}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('linear: n=%.2f, R^2=%.3f, p=%.3g', -a1, {r2_r:.12g}, {p_r:.12g}){log_clause}

set title '{source_label}: RAW (log-log, {normtag}) - log10(sigma) vs log10(S) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Imag Conductivity)'
a2={a_s:.12g}
b2={b_s:.12g}
f2(x)=a2*x+b2
plot '{csv_name}' using (column("logS")):(column("logSigma")) with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('linear: m=%.2f, R^2=%.3f, p=%.3g', a2, {r2_s:.12g}, {p_s:.12g})

unset multiplot
"""
    (results_dir / "plot_raw_loglog.gp").write_text(gp)

def write_plot_real_1panel(results_dir: Path, csv_name: str, out_png: str,
                           title: str, xcol: str, ycol: str,
                           yfit_col: str|None,
                           xlabel: str, ylabel: str,
                           legend_fit: str):
    """
    Generic real-space (no logs) single-panel plot:
      y vs x with optional fitted curve column yfit_col.
    """
    fit_part = ""
    if yfit_col:
        fit_part = f", \\\n     '{csv_name}' using (column(\"{xcol}\")):(column(\"{yfit_col}\")) with lines dt 2 lw 4 title '{legend_fit}'"
    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,26'
set output '{out_png}'
set grid
set key top right
set title '{title}'
set xlabel '{xlabel}'
set ylabel '{ylabel}'
plot '{csv_name}' using (column("{xcol}")):(column("{ycol}")) with points pt 7 ps 1.7 title 'data'{fit_part}
"""
    (results_dir / "plot_real_1panel.gp").write_text(gp)

# =============================================================================
# DISCOVERY
# =============================================================================

def discover_run_roots_fast(base: Path):
    roots = set()
    for a in base.rglob(ANALYSIS_DIR):
        if not a.is_dir():
            continue
        run = a.parent
        if (run / RAW_DIR).is_dir():
            roots.add(run)
    return sorted(roots)

# =============================================================================
# EXCEL WRITERS
# =============================================================================

def write_run_excel(results_dir: Path, xlsx_name: str, df_data: pd.DataFrame,
                    df_fit: pd.DataFrame, df_meta: pd.DataFrame):
    out = results_dir / xlsx_name
    try:
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            df_data.to_excel(w, sheet_name="Data", index=False)
            df_fit.to_excel(w, sheet_name="FitSummary", index=False)
            df_meta.to_excel(w, sheet_name="Meta", index=False)
    except Exception as e:
        print(f"[WARN] Failed to write Excel {out}: {e}", flush=True)

def write_global_excel(base: Path, xlsx_name: str, df_summary: pd.DataFrame, df_allfits: pd.DataFrame|None):
    out = base / xlsx_name
    try:
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            df_summary.to_excel(w, sheet_name="Summary", index=False)
            if df_allfits is not None and len(df_allfits) > 0:
                df_allfits.to_excel(w, sheet_name="AllFitSummaries", index=False)
    except Exception as e:
        print(f"[WARN] Failed to write global Excel {out}: {e}", flush=True)

# =============================================================================
# ONE-RUN PROCESSOR
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
                poly_raw: bool,
                normalize: bool,
                do_real_fits: bool,
                do_real_logistic: bool):
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return None, None  # (summary_row, fit_df_for_global)

    source_label = pretty_source_name(run_root.name)
    results_dir.mkdir(parents=True, exist_ok=True)
    if clean_results:
        clean_results_dir(results_dir)

    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df, sheet_filter=sheet_filter)
    tag = str(target_freq).replace(".","p")
    normtag = "norm" if normalize else "unorm"

    if sip_df.empty:
        return ({"run": str(run_root), "status": "FAIL", "reason": "No SIP sheets readable/target freq missing", "points": 0,
                 "source": source_label, "normalize": normalize}, None)

    maps, used_sheet, overlap_cnt, block_rows = read_logbook_bestblock_by_overlap(
        log_xlsx, sip_df["measurement_norm"].tolist()
    )
    S_map = maps["S"]
    phi_map = maps["phi"]

    # map saturation and porosity
    sip_df["S"] = sip_df["measurement_norm"].map(S_map).astype(float)
    sip_df["phi"] = sip_df["measurement_norm"].map(phi_map).astype(float)

    # theta (volumetric moisture) if phi available
    sip_df["theta"] = sip_df.apply(
        lambda r: float(r["S"] * r["phi"]) if (np.isfinite(r["S"]) and np.isfinite(r["phi"])) else np.nan,
        axis=1
    )

    # logs (for log-log fits and stored outputs)
    sip_df["logS"]     = sip_df["S"].apply(safe_log10)
    sip_df["logRho"]   = sip_df["rho"].apply(safe_log10)
    sip_df["logSigma"] = sip_df["sigma_imag"].apply(safe_log10)
    sip_df["logTheta"] = sip_df["theta"].apply(safe_log10)

    # debug before filtering
    sip_df.copy().to_csv(results_dir / f"debug_before_filter_{tag}Hz_{normtag}.csv", index=False)

    # primary validity for log-log plots: need finite logS, logRho, logSigma
    mask = np.isfinite(sip_df["logS"]) & np.isfinite(sip_df["logRho"]) & np.isfinite(sip_df["logSigma"])
    dropped = sip_df.loc[~mask, ["measurement", "S", "rho", "sigma_imag", "phi", "theta"]].copy()
    if len(dropped) > 0:
        dropped.to_csv(results_dir / f"dropped_rows_{tag}Hz_{normtag}.csv", index=False)
    sip_df = sip_df.loc[mask].copy()

    if len(sip_df) < min_points:
        return ({"run": str(run_root), "status": "FAIL",
                 "reason": f"Not enough matched points (need >= {min_points})",
                 "points": int(len(sip_df)),
                 "logbook_sheet": used_sheet,
                 "logbook_overlap": overlap_cnt,
                 "logbook_block_rows": block_rows,
                 "source": source_label, "normalize": normalize}, None)

    # choose reference point (max S) if normalization is ON
    iref = int(np.nanargmax(sip_df["S"].to_numpy(float)))
    rho_ref = float(sip_df.iloc[iref]["rho"])
    sig_ref = float(sip_df.iloc[iref]["sigma_imag"])
    logRho_ref = safe_log10(rho_ref)
    logSig_ref = safe_log10(sig_ref)

    # =========================
    # BUILD ALL DERIVED COLUMNS
    # =========================

    # RAW log-log uses rho and sigma directly (always valid)
    # For consistency in outputs, define these "raw" columns explicitly:
    sip_df["x_loglog"] = sip_df["logS"]
    sip_df["y_raw_logrho"] = sip_df["logRho"]
    sip_df["y_raw_logsig"] = sip_df["logSigma"]

    # Normalized (classic) indices and VN-style "target" definitions
    if normalize:
        # GOAL/INDEX definitions in log domain:
        sip_df["rho_idx"]   = sip_df["rho"] / rho_ref
        sip_df["logRhoIdx"] = sip_df["logRho"] - logRho_ref

        if ici_mode == "sig_over_sigref":
            # ICI = sig/sig_ref => logICI = logSig - logSig_ref
            sip_df["imag_idx"]   = sip_df["sigma_imag"] / sig_ref
            sip_df["logImagIdx"] = sip_df["logSigma"] - logSig_ref
        else:
            # ICI = sig_ref/sig => logICI = logSig_ref - logSig
            sip_df["imag_idx"]   = sig_ref / sip_df["sigma_imag"]
            sip_df["logImagIdx"] = logSig_ref - sip_df["logSigma"]

        # TARGET (VN style) ALWAYS uses ICI = sig_ref/sig to force downtrend
        sip_df["logRI_vn"]  = sip_df["logRhoIdx"]
        sip_df["logICI_vn"] = logSig_ref - sip_df["logSigma"]

        # Also keep real-space versions of indices
        sip_df["RI_vn_real"]  = sip_df["rho_idx"]  # rho/rho_ref
        sip_df["ICI_vn_real"] = sig_ref / sip_df["sigma_imag"]  # sig_ref/sig
    else:
        # Unnormalized alternative "index/target-like" quantities WITHOUT refs:
        #   RI_unorm  = rho            -> logRI_unorm = logRho
        #   ICI_unorm = 1/sigma_imag   -> logICI_unorm = -logSigma
        # These produce an "ICI-like" downtrend if sigma rises with moisture.
        sip_df["rho_idx"]   = np.nan
        sip_df["logRhoIdx"] = np.nan
        sip_df["imag_idx"]  = np.nan
        sip_df["logImagIdx"]= np.nan

        sip_df["logRI_vn"]  = sip_df["logRho"]
        sip_df["logICI_vn"] = -sip_df["logSigma"]

        sip_df["RI_vn_real"]  = sip_df["rho"]
        sip_df["ICI_vn_real"] = 1.0 / sip_df["sigma_imag"]

    # ======================================
    # FITS: LOG-LOG (KEEP ALL OLD APPROACHES)
    # ======================================

    # RAW log-log linear fits
    slope_r, intercept_r, r2_r, p_r, ss_res_lin_r, aic_lin_r = linear_fit(sip_df["logS"], sip_df["logRho"])
    slope_s, intercept_s, r2_s, p_s, ss_res_lin_s, aic_lin_s = linear_fit(sip_df["logS"], sip_df["logSigma"])
    lin_rho = {"slope": slope_r, "intercept": intercept_r, "r2": r2_r, "p": p_r, "aic": aic_lin_r}
    lin_sig = {"slope": slope_s, "intercept": intercept_s, "r2": r2_s, "p": p_s, "aic": aic_lin_s}

    # "INDEX/GOAL" through-origin only meaningful in normalized mode; in unnormalized we compute through-origin
    # on our alternative logRI_vn/logICI_vn just to keep the plot available (and label it clearly).
    if normalize:
        a_ri,  r2_ri,  _ = linear_fit_through_origin(sip_df["logS"], sip_df["logRhoIdx"])
        a_ici, r2_ici, _ = linear_fit_through_origin(sip_df["logS"], sip_df["logImagIdx"])
        n_idx = -a_ri
        if ici_mode == "sig_over_sigref":
            p_idx = +a_ici
        else:
            p_idx = -a_ici
        sip_df["logRhoIdx_fit_to"] = a_ri  * sip_df["logS"]
        sip_df["logImagIdx_fit_to"]= a_ici * sip_df["logS"]
    else:
        # unnormalized analog:
        a_ri,  r2_ri,  _ = linear_fit_through_origin(sip_df["logS"], sip_df["logRI_vn"])
        a_ici, r2_ici, _ = linear_fit_through_origin(sip_df["logS"], sip_df["logICI_vn"])
        n_idx = -a_ri
        p_idx = -a_ici
        sip_df["logRhoIdx_fit_to"] = np.nan
        sip_df["logImagIdx_fit_to"]= np.nan

    # TARGET (VN or unnormalized analog) with intercept
    slope_RI_vn,  int_RI_vn,  r2_RI_vn,  p_RI_vn,  _, _ = linear_fit(sip_df["logS"], sip_df["logRI_vn"])
    slope_ICI_vn, int_ICI_vn, r2_ICI_vn, p_ICI_vn, _, _ = linear_fit(sip_df["logS"], sip_df["logICI_vn"])
    n_vn = -slope_RI_vn
    p_vn = -slope_ICI_vn

    sip_df["logRI_vn_fit"]  = slope_RI_vn  * sip_df["logS"] + int_RI_vn
    sip_df["logICI_vn_fit"] = slope_ICI_vn * sip_df["logS"] + int_ICI_vn

    # Optional VN-style subset S>=vn_smin (applies equally to normalized and unnormalized analog)
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

    # Polynomial fits on TARGET curves (log-log)
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

    # Logistic fit on RAW logRho vs logS (optional, as before)
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

    # Diagnostics curves and residuals in log-log domain
    logS_arr = sip_df["logS"].to_numpy(float)
    logR_arr = sip_df["logRho"].to_numpy(float)
    lin_pred = slope_r*logS_arr + intercept_r
    sip_df["resid_lin_rho"] = logR_arr - lin_pred
    sip_df["logRho_linfit"] = lin_pred

    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        log_pred = logistic4(logS_arr, L, k, x0, c)
        sip_df["resid_log_rho"] = logR_arr - log_pred
        sip_df["logRho_logfit"] = log_pred
    else:
        sip_df["resid_log_rho"] = np.nan
        sip_df["logRho_logfit"] = np.nan

    # ===========================================
    # REAL-SPACE FITS: rho vs S, sigma vs S/theta
    # ===========================================

    # Default to no real fits unless enabled
    realfits = {}

    sip_df["rho_fit_real_lin_S"] = np.nan
    sip_df["sig_fit_real_lin_S"] = np.nan
    sip_df["sig_fit_real_lin_theta"] = np.nan

    sip_df["rho_fit_real_poly_S"] = np.nan
    sip_df["sig_fit_real_poly_S"] = np.nan
    sip_df["sig_fit_real_poly_theta"] = np.nan

    sip_df["rho_fit_real_logistic_S"] = np.nan
    sip_df["rho_resid_real_lin_S"] = np.nan
    sip_df["sig_resid_real_lin_S"] = np.nan

    # For real fits we require finite S and finite y (rho/sigma)
    if do_real_fits:
        xS = sip_df["S"].to_numpy(float)
        yR = sip_df["rho"].to_numpy(float)
        yG = sip_df["sigma_imag"].to_numpy(float)

        # linear real: rho ~ a*S + b
        ar, br, r2r_real, pr_real, ssr_real, aicr_real = linear_fit(xS, yR)
        sip_df["rho_fit_real_lin_S"] = ar*sip_df["S"] + br
        sip_df["rho_resid_real_lin_S"] = sip_df["rho"] - sip_df["rho_fit_real_lin_S"]

        # linear real: sigma ~ a*S + b
        ag, bg, r2g_real, pg_real, ssg_real, aicg_real = linear_fit(xS, yG)
        sip_df["sig_fit_real_lin_S"] = ag*sip_df["S"] + bg
        sip_df["sig_resid_real_lin_S"] = sip_df["sigma_imag"] - sip_df["sig_fit_real_lin_S"]

        realfits["rho_lin_S"] = {"a": ar, "b": br, "R2": r2r_real, "p": pr_real, "AIC": aicr_real}
        realfits["sig_lin_S"] = {"a": ag, "b": bg, "R2": r2g_real, "p": pg_real, "AIC": aicg_real}

        # poly real if requested (re-use --poly-deg)
        if poly_deg and poly_deg >= 2:
            cRp, _, r2Rp, _ = poly_fit_predict(xS, yR, poly_deg)
            cGp, _, r2Gp, _ = poly_fit_predict(xS, yG, poly_deg)
            if cRp is not None:
                sip_df["rho_fit_real_poly_S"] = np.polyval(cRp, sip_df["S"].to_numpy(float))
            if cGp is not None:
                sip_df["sig_fit_real_poly_S"] = np.polyval(cGp, sip_df["S"].to_numpy(float))
            realfits["rho_poly_S"] = {"deg": poly_deg, "R2": r2Rp}
            realfits["sig_poly_S"] = {"deg": poly_deg, "R2": r2Gp}

        # sigma vs theta if theta exists
        mth = np.isfinite(sip_df["theta"].to_numpy(float))
        if int(np.sum(mth)) >= min_points:
            xt = sip_df.loc[mth, "theta"].to_numpy(float)
            yg = sip_df.loc[mth, "sigma_imag"].to_numpy(float)
            at, bt, r2t, pt, sst, aict = linear_fit(xt, yg)
            sip_df.loc[mth, "sig_fit_real_lin_theta"] = at*sip_df.loc[mth, "theta"] + bt
            realfits["sig_lin_theta"] = {"a": at, "b": bt, "R2": r2t, "p": pt, "AIC": aict}

            if poly_deg and poly_deg >= 2 and int(np.sum(mth)) >= (poly_deg + 1):
                cT, _, r2Tp, _ = poly_fit_predict(xt, yg, poly_deg)
                if cT is not None:
                    sip_df.loc[mth, "sig_fit_real_poly_theta"] = np.polyval(cT, sip_df.loc[mth, "theta"].to_numpy(float))
                realfits["sig_poly_theta"] = {"deg": poly_deg, "R2": r2Tp}

        # optional real-space logistic on rho vs S
        if do_real_logistic and len(sip_df) >= logistic_min_points:
            lf_real, st_real = logistic_fit_noscipy(
                sip_df["S"], sip_df["rho"],
                fast=logistic_fast,
                max_seconds=(logistic_max_seconds if logistic_max_seconds and logistic_max_seconds > 0 else 0.0)
            )
            if lf_real is not None:
                okr, rr = accept_logistic(r2r_real, aicr_real, lf_real)
                if okr:
                    sip_df["rho_fit_real_logistic_S"] = logistic4(
                        sip_df["S"].to_numpy(float), lf_real["L"], lf_real["k"], lf_real["x0"], lf_real["c"]
                    )
                    realfits["rho_logistic_S"] = {"status": st_real, "accepted": True, "reason": rr, **lf_real}
                else:
                    realfits["rho_logistic_S"] = {"status": st_real, "accepted": False, "reason": rr, **lf_real}

    # ======================
    # WRITE JOINED CSV
    # ======================

    csv_path = results_dir / f"joined_{tag}Hz_{normtag}.csv"
    cols = [
        "measurement","measurement_norm","freq_hz",
        "S","phi","theta",
        "rho","sigma_imag",
        "logS","logTheta",
        "logRho","logSigma",

        # normalized indices if available, else NaN
        "rho_idx","imag_idx","logRhoIdx","logImagIdx",

        # target/vn columns (normalized or unnormalized analog)
        "logRI_vn","logICI_vn",
        "logRI_vn_fit","logICI_vn_fit",
        "logRI_vn_fit_smin","logICI_vn_fit_smin",
        "logRI_vn_polyfit","logICI_vn_polyfit",
        "logRI_vn_polyfit_smin","logICI_vn_polyfit_smin",

        # raw poly
        "logRho_polyfit",

        # log-log diagnostics
        "logRho_linfit","logRho_logfit",
        "resid_lin_rho","resid_log_rho",

        # real fits
        "rho_fit_real_lin_S","sig_fit_real_lin_S",
        "rho_fit_real_poly_S","sig_fit_real_poly_S",
        "rho_fit_real_logistic_S",
        "sig_fit_real_lin_theta","sig_fit_real_poly_theta",
        "rho_resid_real_lin_S","sig_resid_real_lin_S",
    ]
    keep = [c for c in cols if c in sip_df.columns]
    sip_df[keep].to_csv(csv_path, index=False)

    # ======================
    # REPORT (TXT)
    # ======================

    report = results_dir / f"fit_report_{tag}Hz_{normtag}.txt"
    report_lines = []
    report_lines += [
        f"Run root: {run_root}",
        f"Source (folder): {source_label}",
        f"Normalize: {normalize}",
        f"SIP file: {sip_xlsx}",
        f"Logbook file: {log_xlsx}",
        f"Logbook sheet used: {used_sheet}",
        f"Logbook block rows (usable): {block_rows}",
        f"Overlap with SIP measurement IDs: {overlap_cnt}",
        "",
        f"Matched points used: {len(sip_df)}",
        f"Reference point (max S): {sip_df.iloc[iref]['measurement']}  S_ref={sip_df.iloc[iref]['S']:.6g}",
        f"rho_ref={rho_ref:.6g}, sigma_ref={sig_ref:.6g}",
        "",
        f"GOAL/INDEX ici_mode: {ici_mode}",
        "",
        "LOG-LOG FITS (x=logS):",
        f"  RAW logRho=a*logS+b: a={slope_r:.6g}, b={intercept_r:.6g}, n(raw)={-slope_r:.6g}, R2={r2_r:.6g}, p={p_r:.3g}",
        f"  RAW logSigma=a*logS+b: a={slope_s:.6g}, b={intercept_s:.6g}, R2={r2_s:.6g}, p={p_s:.3g}",
        "",
        f"TARGET (VN-like; WITH intercept) exponents:",
        f"  n_vn={n_vn:.6g} (R2={r2_RI_vn:.6g}, p={p_RI_vn:.3g})",
        f"  p_vn={p_vn:.6g} (R2={r2_ICI_vn:.6g}, p={p_ICI_vn:.3g})",
        "",
        "INDEX/GOAL (through-origin):",
        f"  n(index)={n_idx:.6g} (R2={r2_ri:.6g})",
        f"  p(index)={p_idx:.6g} (R2={r2_ici:.6g})",
        "",
        f"Logistic (log-log rho): status={log_status}, accepted={bool(log_fit_acc is not None)}, reason={log_accept_reason}",
        "",
    ]
    if vn_smin and vn_smin > 0:
        report_lines += [
            f"EXTRA TARGET subset: S >= {vn_smin:g}",
            f"  n_vn_smin={n_vn_smin:.6g} (R2={r2_RI_vn_smin:.6g})",
            f"  p_vn_smin={p_vn_smin:.6g} (R2={r2_ICI_vn_smin:.6g})",
            "",
        ]
    if poly_deg and poly_deg >= 2:
        report_lines += [
            f"POLY (log-log) TARGET deg {poly_deg}:",
            f"  R2_RI_vn_poly={r2_RI_vn_poly:.6g}",
            f"  R2_ICI_vn_poly={r2_ICI_vn_poly:.6g}",
            "",
        ]
        if vn_smin and vn_smin > 0:
            report_lines += [
                f"POLY (log-log) TARGET deg {poly_deg}, S >= {vn_smin:g}:",
                f"  R2_RI_vn_poly_smin={r2_RI_vn_poly_smin:.6g}",
                f"  R2_ICI_vn_poly_smin={r2_ICI_vn_poly_smin:.6g}",
                "",
            ]
        if poly_raw:
            report_lines += [f"POLY (log-log) RAW logRho deg {poly_deg}: R2_raw_poly={r2_raw_poly:.6g}", ""]

    if do_real_fits and realfits:
        report_lines += ["REAL-SPACE FITS (no logs):"]
        for k, v in realfits.items():
            report_lines.append(f"  {k}: {v}")
        report_lines.append("")

    report.write_text("\n".join(report_lines))

    # ======================
    # EXCEL OUTPUT (PER RUN)
    # ======================

    df_meta = pd.DataFrame([{
        "run_root": str(run_root),
        "source": source_label,
        "normalize": normalize,
        "normtag": normtag,
        "target_freq": target_freq,
        "sip_xlsx": str(sip_xlsx),
        "log_xlsx": str(log_xlsx),
        "logbook_sheet": used_sheet,
        "logbook_overlap": overlap_cnt,
        "logbook_block_rows": block_rows,
        "reference_measurement": str(sip_df.iloc[iref]["measurement"]),
        "S_ref": float(sip_df.iloc[iref]["S"]),
        "rho_ref": rho_ref,
        "sigma_ref": sig_ref
    }])

    # Fit summary sheet: one row with many columns (easy to aggregate globally)
    fit_row = {
        "run_root": str(run_root),
        "source": source_label,
        "normalize": normalize,
        "normtag": normtag,
        "freq_hz": target_freq,
        "points": int(len(sip_df)),
        "ici_mode": ici_mode,

        # RAW log-log
        "raw_loglog_slope_logRho_vs_logS": slope_r,
        "raw_loglog_int_logRho_vs_logS": intercept_r,
        "raw_loglog_R2_logRho": r2_r,
        "raw_loglog_p_logRho": p_r,

        "raw_loglog_slope_logSigma_vs_logS": slope_s,
        "raw_loglog_int_logSigma_vs_logS": intercept_s,
        "raw_loglog_R2_logSigma": r2_s,
        "raw_loglog_p_logSigma": p_s,

        # TARGET log-log
        "target_loglog_slope_logRI": slope_RI_vn,
        "target_loglog_int_logRI": int_RI_vn,
        "target_loglog_R2_logRI": r2_RI_vn,
        "target_loglog_n": n_vn,

        "target_loglog_slope_logICI": slope_ICI_vn,
        "target_loglog_int_logICI": int_ICI_vn,
        "target_loglog_R2_logICI": r2_ICI_vn,
        "target_loglog_p": p_vn,

        # INDEX
        "index_through_origin_n": n_idx,
        "index_through_origin_R2_RI": r2_ri,
        "index_through_origin_p": p_idx,
        "index_through_origin_R2_ICI": r2_ici,

        # Logistic log-log on rho
        "logistic_loglog_attempted": bool(log_fit is not None),
        "logistic_loglog_status": log_status,
        "logistic_loglog_accepted": bool(log_fit_acc is not None),
        "logistic_loglog_reason": log_accept_reason,
        "logistic_loglog_r2": (log_fit_acc.get("r2") if log_fit_acc else np.nan),
        "logistic_loglog_p_model": (log_fit_acc.get("p_model") if log_fit_acc else np.nan),
        "logistic_loglog_aic": (log_fit_acc.get("aic") if log_fit_acc else np.nan),
    }

    # poly metrics
    if poly_deg and poly_deg >= 2:
        fit_row.update({
            "poly_deg": poly_deg,
            "target_poly_R2_RI": r2_RI_vn_poly,
            "target_poly_R2_ICI": r2_ICI_vn_poly,
            "target_poly_R2_RI_Smin": r2_RI_vn_poly_smin,
            "target_poly_R2_ICI_Smin": r2_ICI_vn_poly_smin,
            "raw_poly_R2_logRho": r2_raw_poly if poly_raw else np.nan
        })
    else:
        fit_row.update({"poly_deg": 0})

    # real metrics if present
    if do_real_fits:
        # store in a flat way (safe if missing)
        rl = realfits.get("rho_lin_S", {})
        sl = realfits.get("sig_lin_S", {})
        st = realfits.get("sig_lin_theta", {})
        fit_row.update({
            "real_rho_lin_S_a": rl.get("a", np.nan),
            "real_rho_lin_S_b": rl.get("b", np.nan),
            "real_rho_lin_S_R2": rl.get("R2", np.nan),
            "real_sig_lin_S_a": sl.get("a", np.nan),
            "real_sig_lin_S_b": sl.get("b", np.nan),
            "real_sig_lin_S_R2": sl.get("R2", np.nan),
            "real_sig_lin_theta_a": st.get("a", np.nan),
            "real_sig_lin_theta_b": st.get("b", np.nan),
            "real_sig_lin_theta_R2": st.get("R2", np.nan),
        })
        if poly_deg and poly_deg >= 2:
            rp = realfits.get("rho_poly_S", {})
            sp = realfits.get("sig_poly_S", {})
            tp = realfits.get("sig_poly_theta", {})
            fit_row.update({
                "real_rho_poly_S_R2": rp.get("R2", np.nan),
                "real_sig_poly_S_R2": sp.get("R2", np.nan),
                "real_sig_poly_theta_R2": tp.get("R2", np.nan),
            })
        lg = realfits.get("rho_logistic_S", {})
        if lg:
            fit_row.update({
                "real_rho_logistic_S_status": lg.get("status", ""),
                "real_rho_logistic_S_accepted": lg.get("accepted", False),
                "real_rho_logistic_S_reason": lg.get("reason", ""),
                "real_rho_logistic_S_R2": lg.get("r2", np.nan),
                "real_rho_logistic_S_aic": lg.get("aic", np.nan),
                "real_rho_logistic_S_p_model": lg.get("p_model", np.nan),
                "real_rho_logistic_S_L": lg.get("L", np.nan),
                "real_rho_logistic_S_k": lg.get("k", np.nan),
                "real_rho_logistic_S_x0": lg.get("x0", np.nan),
                "real_rho_logistic_S_c": lg.get("c", np.nan),
            })

    df_fit = pd.DataFrame([fit_row])

    # Write run excel
    write_run_excel(
        results_dir,
        xlsx_name=f"all_values_{tag}Hz_{normtag}.xlsx",
        df_data=sip_df.copy(),
        df_fit=df_fit,
        df_meta=df_meta
    )

    # ======================
    # PLOTS (LOG-LOG, KEEP OLD + ADD NORM TAG IN NAME)
    # ======================

    # RAW log-log
    out_raw = f"exponents_raw_{tag}Hz_{normtag}.png"
    write_plot_raw_loglog(results_dir, csv_path.name, out_raw, target_freq, source_label, lin_rho, lin_sig, log_fit_acc, normtag)

    # TARGET log-log (VN or unnormalized analog)
    out_tgt = f"target_exponents_{tag}Hz_{normtag}.png"
    # reuse a generic 2-panel writer
    write_plot_target_loglog(
        results_dir, csv_path.name, out_tgt, target_freq, source_label,
        n_val=n_vn, r2a=r2_RI_vn, p_val=p_vn, r2b=r2_ICI_vn,
        col_y1="logRI_vn", col_y1_fit="logRI_vn_fit",
        col_y2="logICI_vn", col_y2_fit="logICI_vn_fit",
        title1="TARGET (log-log)", title2="TARGET (log-log)",
        ylab1="log10(RI-like)", ylab2="log10(ICI-like)",
        legend1=f"{label_fmt(target_freq)}Hz: n = %.2f, R^2 = %.3f",
        legend2=f"{label_fmt(target_freq)}Hz: p = %.2f, R^2 = %.3f",
    )

    # optional poly target plots (log-log)
    if poly_deg and poly_deg >= 2:
        # 2-panel plot overlay with poly columns is easy to create by cloning gp,
        # but to keep code compact we generate dedicated plots via 1-panel scripts per y.
        # Instead, we keep your earlier naming convention and just plot poly curves as additional file:
        out_tgt_poly = f"target_exponents_polydeg{poly_deg}_{tag}Hz_{normtag}.png"
        gp = f"""reset
set datafile separator ','
set datafile missing ''
{_gp_term()}set output '{out_tgt_poly}'
set multiplot layout 2,1
set grid
set key top right
set xlabel 'log10(Saturation)'

set title '{source_label}: TARGET poly(deg {poly_deg}) (log-log, {normtag}) - RI-like at {target_freq:g} Hz'
set ylabel 'log10(RI-like)'
plot '{csv_path.name}' using (column("logS")):(column("logRI_vn")) with points pt 7 ps 1.7 title 'data', \\
     '{csv_path.name}' using (column("logS")):(column("logRI_vn_polyfit")) with lines dt 1 lw 4 title sprintf('poly deg {poly_deg}: R^2=%.3f', {r2_RI_vn_poly if np.isfinite(r2_RI_vn_poly) else float("nan")})

set title '{source_label}: TARGET poly(deg {poly_deg}) (log-log, {normtag}) - ICI-like at {target_freq:g} Hz'
set ylabel 'log10(ICI-like)'
plot '{csv_path.name}' using (column("logS")):(column("logICI_vn")) with points pt 7 ps 1.7 title 'data', \\
     '{csv_path.name}' using (column("logS")):(column("logICI_vn_polyfit")) with lines dt 1 lw 4 title sprintf('poly deg {poly_deg}: R^2=%.3f', {r2_ICI_vn_poly if np.isfinite(r2_ICI_vn_poly) else float("nan")})

unset multiplot
"""
        (results_dir / "plot_target_poly_loglog.gp").write_text(gp)

    # optional RAW poly logRho plot (log-log)
    if poly_raw and (poly_deg and poly_deg >= 2):
        out_raw_poly = f"exponents_raw_polydeg{poly_deg}_{tag}Hz_{normtag}.png"
        gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,26'
set output '{out_raw_poly}'
set grid
set key top right
set title '{source_label}: RAW poly(deg {poly_deg}) (log-log, {normtag}) - logRho vs logS at {target_freq:g} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
plot '{csv_path.name}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.7 title 'data', \\
     '{csv_path.name}' using (column("logS")):(column("logRho_polyfit")) with lines dt 1 lw 4 title sprintf('poly deg {poly_deg}: R^2=%.3f', {r2_raw_poly if np.isfinite(r2_raw_poly) else float("nan")})
"""
        (results_dir / "plot_raw_poly_loglog.gp").write_text(gp)

    # ======================
    # REAL-SPACE PLOTS
    # ======================
    if do_real_fits:
        # rho vs S (linear)
        write_plot_real_1panel(
            results_dir, csv_path.name,
            out_png=f"real_rho_vs_S_linear_{tag}Hz_{normtag}.png",
            title=f"{source_label}: rho vs S (real, {normtag}) at {target_freq:g} Hz",
            xcol="S", ycol="rho", yfit_col="rho_fit_real_lin_S",
            xlabel="Saturation (S)", ylabel="Resistivity (rho)",
            legend_fit="linear fit"
        )
        # sigma vs S (linear) -> conductivity vs moisture (core goal)
        write_plot_real_1panel(
            results_dir, csv_path.name,
            out_png=f"real_sigma_vs_S_linear_{tag}Hz_{normtag}.png",
            title=f"{source_label}: sigma_imag vs S (real, {normtag}) at {target_freq:g} Hz",
            xcol="S", ycol="sigma_imag", yfit_col="sig_fit_real_lin_S",
            xlabel="Saturation (S)", ylabel="Imag Conductivity (sigma_imag)",
            legend_fit="linear fit"
        )
        # sigma vs theta (linear) if theta exists
        if np.isfinite(sip_df["theta"].to_numpy(float)).any():
            write_plot_real_1panel(
                results_dir, csv_path.name,
                out_png=f"real_sigma_vs_theta_linear_{tag}Hz_{normtag}.png",
                title=f"{source_label}: sigma_imag vs theta (real, {normtag}) at {target_freq:g} Hz",
                xcol="theta", ycol="sigma_imag", yfit_col="sig_fit_real_lin_theta",
                xlabel="Volumetric moisture (theta = S*phi)", ylabel="Imag Conductivity (sigma_imag)",
                legend_fit="linear fit"
            )
        # poly overlays if enabled
        if poly_deg and poly_deg >= 2:
            if np.isfinite(sip_df["rho_fit_real_poly_S"].to_numpy(float)).any():
                write_plot_real_1panel(
                    results_dir, csv_path.name,
                    out_png=f"real_rho_vs_S_polydeg{poly_deg}_{tag}Hz_{normtag}.png",
                    title=f"{source_label}: rho vs S (real, poly deg {poly_deg}, {normtag}) at {target_freq:g} Hz",
                    xcol="S", ycol="rho", yfit_col="rho_fit_real_poly_S",
                    xlabel="Saturation (S)", ylabel="Resistivity (rho)",
                    legend_fit=f"poly deg {poly_deg}"
                )
            if np.isfinite(sip_df["sig_fit_real_poly_S"].to_numpy(float)).any():
                write_plot_real_1panel(
                    results_dir, csv_path.name,
                    out_png=f"real_sigma_vs_S_polydeg{poly_deg}_{tag}Hz_{normtag}.png",
                    title=f"{source_label}: sigma_imag vs S (real, poly deg {poly_deg}, {normtag}) at {target_freq:g} Hz",
                    xcol="S", ycol="sigma_imag", yfit_col="sig_fit_real_poly_S",
                    xlabel="Saturation (S)", ylabel="Imag Conductivity (sigma_imag)",
                    legend_fit=f"poly deg {poly_deg}"
                )
            if np.isfinite(sip_df["sig_fit_real_poly_theta"].to_numpy(float)).any():
                write_plot_real_1panel(
                    results_dir, csv_path.name,
                    out_png=f"real_sigma_vs_theta_polydeg{poly_deg}_{tag}Hz_{normtag}.png",
                    title=f"{source_label}: sigma_imag vs theta (real, poly deg {poly_deg}, {normtag}) at {target_freq:g} Hz",
                    xcol="theta", ycol="sigma_imag", yfit_col="sig_fit_real_poly_theta",
                    xlabel="Volumetric moisture (theta = S*phi)", ylabel="Imag Conductivity (sigma_imag)",
                    legend_fit=f"poly deg {poly_deg}"
                )
        # logistic rho vs S
        if do_real_logistic and np.isfinite(sip_df["rho_fit_real_logistic_S"].to_numpy(float)).any():
            write_plot_real_1panel(
                results_dir, csv_path.name,
                out_png=f"real_rho_vs_S_logistic_{tag}Hz_{normtag}.png",
                title=f"{source_label}: rho vs S (real, logistic, {normtag}) at {target_freq:g} Hz",
                xcol="S", ycol="rho", yfit_col="rho_fit_real_logistic_S",
                xlabel="Saturation (S)", ylabel="Resistivity (rho)",
                legend_fit="logistic fit"
            )

    # ======================
    # RUN GNUPLOT
    # ======================
    if not dry_run:
        _run_gnuplot("plot_raw_loglog.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_target_loglog.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)

        if poly_deg and poly_deg >= 2:
            _run_gnuplot("plot_target_poly_loglog.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        if poly_raw and (poly_deg and poly_deg >= 2):
            _run_gnuplot("plot_raw_poly_loglog.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)

        if do_real_fits:
            _run_gnuplot("plot_real_1panel.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)
            # NOTE: We overwrite plot_real_1panel.gp multiple times; gnuplot runs only last one if we run once.
            # To avoid that, we run gnuplot immediately after each write when not dry_run.
            # For simplicity, we do a second pass here by directly running gnuplot on each plot by regenerating.
            # Instead of complexity, we just re-run gnuplot by listing all .gp files in Results.
            for gpfile in sorted(results_dir.glob("plot_*.gp")):
                _run_gnuplot(gpfile.name, cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

    # ======================
    # RETURN SUMMARY + FIT DF
    # ======================
    summary = {
        "run": str(run_root),
        "source": source_label,
        "status": "OK",
        "points": int(len(sip_df)),
        "normalize": normalize,
        "normtag": normtag,
        "logbook_sheet": used_sheet,
        "logbook_overlap": overlap_cnt,
        "logbook_block_rows": block_rows,
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
        "poly_raw": bool(poly_raw),
        "R2_raw_poly": float(r2_raw_poly) if np.isfinite(r2_raw_poly) else r2_raw_poly,

        "real_fits": bool(do_real_fits),
        "real_logistic": bool(do_real_logistic),
    }

    return summary, df_fit


def label_fmt(freq: float) -> str:
    # helper for labels inside gnuplot sprintf, avoid too many decimals
    return "{:g}".format(freq)

# =============================================================================
# MAIN CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("top", help="Top folder containing run subfolders")
    ap.add_argument("--freq", type=float, default=0.01, help="Target frequency in Hz")
    ap.add_argument("--sample", type=str, default=None, help="Substring filter on SIP sheet names (debug)")
    ap.add_argument("--dry-run", action="store_true", help="Write scripts/CSVs but do not run gnuplot")
    ap.add_argument("--min-points", type=int, default=4, help="Minimum matched points to accept a run")
    ap.add_argument("--max-df", type=float, default=0.0, help="Max allowed |f_found - f_target| (0 disables)")
    ap.add_argument("--no-logistic", action="store_true", help="Disable logistic fit/plots (log-log rho)")
    ap.add_argument("--logistic-min-points", type=int, default=6, help="Min points needed to attempt logistic")
    ap.add_argument("--logistic-fast", action="store_true", help="Use SciPy-free fast logistic search")
    ap.add_argument("--logistic-max-seconds", type=float, default=8.0, help="Time cap for logistic fit per run")
    ap.add_argument("--gnuplot-timeout", type=int, default=30, help="Timeout seconds per gnuplot call")
    ap.add_argument("--clean-results", action="store_true", help="Delete files in each run's Results/ before writing")

    ap.add_argument("--ici-mode", type=str, default="sigref_over_sig",
                    choices=["sigref_over_sig","sig_over_sigref"],
                    help=("GOAL/INDEX ICI definition in NORMALIZED mode. "
                          "In UNNORMALIZED mode we do not use a ref; we use ICI_unorm = 1/sigma_imag."))

    ap.add_argument("--vn-smin", type=float, default=0.0,
                    help="Extra TARGET fit/plot using only points with S >= vn_smin.")

    ap.add_argument("--poly-deg", type=int, default=0,
                    help="Polynomial degree for extra fits/plots (>=2). 0 disables.")
    ap.add_argument("--poly-raw", action="store_true",
                    help="Also plot polynomial fit for RAW logRho vs logS (requires --poly-deg>=2).")

    # NEW: normalization switch
    ap.add_argument("--normalize", dest="normalize", action="store_true", help="Enable normalization (default).")
    ap.add_argument("--no-normalize", dest="normalize", action="store_false", help="Disable ref-based normalization.")
    ap.set_defaults(normalize=True)

    # NEW: real-space fits/plots switches
    ap.add_argument("--real-fits", action="store_true", help="Enable real-space fits/plots (rho vs S, sigma vs S/theta).")
    ap.add_argument("--real-logistic", action="store_true", help="Also attempt real-space logistic rho vs S (requires --real-fits).")

    args = ap.parse_args()

    base = Path(args.top).expanduser().resolve()
    if not base.exists():
        print(f"Not found: {base}")
        raise SystemExit(2)

    roots = discover_run_roots_fast(base)
    print(f"Found {len(roots)} run folder(s).", flush=True)

    rows = []
    allfits = []
    t0 = time.time()

    for i, run in enumerate(roots, 1):
        print(f"[{i}/{len(roots)}] START {run}", flush=True)
        ti = time.time()
        try:
            summary, df_fit = process_run(
                run_root=run,
                target_freq=args.freq,
                sheet_filter=args.sample,
                dry_run=args.dry_run,
                min_points=args.min_points,
                max_df=args.max_df,
                do_logistic=(not args.no_logistic),
                logistic_min_points=args.logistic_min_points,
                logistic_fast=args.logistic_fast,
                logistic_max_seconds=args.logistic_max_seconds,
                gnuplot_timeout=args.gnuplot_timeout,
                clean_results=args.clean_results,
                ici_mode=args.ici_mode,
                vn_smin=args.vn_smin,
                poly_deg=args.poly_deg,
                poly_raw=args.poly_raw,
                normalize=args.normalize,
                do_real_fits=args.real_fits,
                do_real_logistic=(args.real_logistic and args.real_fits),
            )
            if summary is None:
                continue
            rows.append(summary)
            if df_fit is not None:
                allfits.append(df_fit)

            dt = time.time() - ti
            msg = (f"[{i}/{len(roots)}] OK points={summary.get('points',0)} "
                   f"norm={summary.get('normalize',True)} "
                   f"n_vn={summary.get('n_vn',np.nan):.3g} p_vn={summary.get('p_vn',np.nan):.3g} "
                   f"({dt:.1f}s)")
            print(msg, flush=True)

        except subprocess.TimeoutExpired:
            rows.append({"run": str(run), "status": "FAIL", "reason": "gnuplot timeout", "points": 0, "normalize": args.normalize})
            print(f"[{i}/{len(roots)}] FAIL gnuplot timeout", flush=True)
        except Exception as e:
            rows.append({"run": str(run), "status": "FAIL", "reason": str(e), "points": 0, "normalize": args.normalize})
            print(f"[{i}/{len(roots)}] FAIL {e}", flush=True)

    tag = str(args.freq).replace(".","p")
    normtag = "norm" if args.normalize else "unorm"

    if rows:
        df_summary = pd.DataFrame(rows)
        summary_csv = base / f"Results_Summary_{tag}Hz_{normtag}.csv"
        df_summary.to_csv(summary_csv, index=False)
        print(f"\nSummary written: {summary_csv}", flush=True)

        df_allfits = pd.concat(allfits, ignore_index=True) if allfits else None
        write_global_excel(
            base,
            xlsx_name=f"Results_All_{tag}Hz_{normtag}.xlsx",
            df_summary=df_summary,
            df_allfits=df_allfits
        )
        print(f"Global Excel written: {base / f'Results_All_{tag}Hz_{normtag}.xlsx'}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)

if __name__ == "__main__":
    main()
