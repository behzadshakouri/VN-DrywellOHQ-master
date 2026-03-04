#!/usr/bin/env python3
"""
auto_sip_pipeline_v27_noscipy.py  (v27++ : FULLY DOCUMENTED / COMMENTED)

WHAT YOU ASKED FOR (IMPLEMENTED)
--------------------------------
Keep *all* previous v27 approaches/plots/logic, and ADD:

1) Run **ALL methods** in **BOTH normalization modes**:
     - NORM  : reference-based indices (rho_ref, sigma_ref; RI/ICI) + VN target.
     - UNORM : no reference normalization (direct rho, sigma_imag).

   You can do this in ONE run using:
     --both-modes
   And you can "clean once, then run both" with:
     --clean-results --both-modes

2) Run **ALL methods** in **BOTH spaces**:
     - LOGLOG: fits performed in log10 space (existing behavior + new additions)
     - REAL  : fits performed in real space (new)

3) Run **ALL fit methods** (where applicable) in each (mode, space):
     - Linear (OLS, with intercept)
     - Polynomial (deg>=2) if --poly-deg >=2
     - 4-parameter logistic (SciPy-free grid/coordinate search) if enabled
       (Your logistic4 is already 4-parameter; now applied broadly.)

4) Output organization: filenames clearly encode:
     mode:  NORM / UNORM
     space: LOGLOG / REAL
     y-x:   e.g., RHO_vs_S, SIG_vs_S, RI_vs_S, ICI_vs_S, etc.
     method: LINEAR / POLYdegD / LOGISTIC4
   so you can find results instantly.

5) Excel outputs (per run folder, AND global):
   - Per run: Results/all_values_<tag>Hz.xlsx
       Sheets:
         - data__NORM__LOGLOG
         - data__NORM__REAL
         - data__UNORM__LOGLOG
         - data__UNORM__REAL
         - fits__NORM
         - fits__UNORM
   - Global: <top>/Results_Summary_<tag>Hz.xlsx and Results_Summary_<tag>Hz.csv

6) "Final goal = relation between electric conductivity and moisture":
   Implemented as systematic fits/plots for:
     sigma_imag  vs S     (REAL + LOGLOG, NORM + UNORM)
   plus optional (if you later add theta column) you can extend similarly.
   (Right now moisture input is saturation S from logbook mapping.)

NOTES / CONVENTIONS
-------------------
- "S" is saturation (fraction 0..1).
- "rho" is resistivity at the target frequency.
- "sigma_imag" is imaginary conductivity at the target frequency (units as in your SIP output).
- "NORM" mode:
    rho_ref = rho at max S
    sig_ref = sigma_imag at max S
    RI = rho / rho_ref
    GOAL/INDEX ICI uses --ici-mode:
       sigref_over_sig (default): ICI = sig_ref / sig
       sig_over_sigref           : ICI = sig / sig_ref
    TARGET "VN style" always uses ICI_vn = sig_ref / sig (downtrend by definition), independent of --ici-mode.
- "UNORM" mode:
    We do NOT compute RI/ICI indices for the *main* unnormalized conductivity relation.
    Instead we fit:
      rho vs S, sigma_imag vs S  (REAL + LOGLOG) with linear/poly/logistic4.
    (You still keep original RAW log-log plots; here they are duplicated with UNORM tags.)

DEPENDENCIES
------------
- numpy, pandas
- gnuplot installed for plots

USAGE
-----
# Clean once, run BOTH modes, all spaces, polynomial, and logistic4 broadly
python3 auto_sip_pipeline_v27_noscipy.py "/path/to/SIP" \
  --freq 0.01 --min-points 4 \
  --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
  --clean-results --both-modes \
  --ici-mode sig_over_sigref \
  --poly-deg 2 --poly-raw \
  --logistic4

TIP
---
If you want EXACTLY "clean then do both norm and unorm" in one invocation:
  use --clean-results --both-modes
"""

import math
import argparse
import subprocess
import time
import re
from pathlib import Path
from dataclasses import dataclass, asdict

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

# Logistic bounds / acceptance (unchanged logic)
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
    # Put S-1 into enhanced text braces for gnuplot: S_{-1}
    s = re.sub(r'\bS-(\d+)\b', r'S_{-\1}', s)
    return s

def tag_freq(freq: float) -> str:
    return str(freq).replace(".", "p")

def mode_suffix(normalize: bool) -> str:
    return "NORM" if normalize else "UNORM"

def space_suffix(is_log: bool) -> str:
    return "LOGLOG" if is_log else "REAL"

def fname(prefix: str, normalize: bool, is_log: bool, tag: str, extra: str = "", ext: str = ".png") -> str:
    """
    Build consistent filenames encoding mode/space and optional extra info.
    Example:
      cond_vs_S__REAL__UNORM__LINEAR+LOGISTIC4_0p01Hz.png
    """
    base = f"{prefix}__{space_suffix(is_log)}__{mode_suffix(normalize)}"
    if extra:
        base += f"__{extra}"
    base += f"_{tag}Hz{ext}"
    return base


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
    """
    OLS y = slope*x + intercept
    Returns: slope, intercept, r2, pvalue, ss_res, aic
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
    Forced through origin y = a*x
    Returns: a, r2, ss_res
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
    Returns (coeffs_high_to_low, yhat_on_cleaned_x, r2, ss_res).
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
# 4-PARAM LOGISTIC FIT (NO SCIPY)  -- ALREADY 4-PARAM, NOW USED BROADLY
# =============================================================================

def logistic4(x, L, k, x0, c):
    """
    4-parameter logistic:
        f(x)=c + L/(1+exp(-k*(x-x0)))
    """
    z = k*(x - x0)
    z = np.clip(z, -60.0, 60.0)
    return c + L/(1.0 + np.exp(-z))

def _logistic_init_guess(x, y):
    """
    Heuristic init:
      c ~ min(y)
      L ~ max(y)-min(y)
      x0 ~ x at y ~ mid
      k ~ slope near mid / L
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
    Fast SciPy-free 4-param logistic via coarse coordinate grid search.
    Returns (fit_dict, status).
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    n = len(x)
    if n < 5:
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
    p_model = f_test_pvalue(ss_res, ss_tot, df_model=3, n=n)  # effective predictors
    aic = aic_gaussian(ss_res, n, k_params=4)

    return {
        "L": float(Lb), "k": float(kb), "x0": float(x0b), "c": float(cb),
        "r2": float(r2),
        "p_model": float(p_model) if np.isfinite(p_model) else p_model,
        "aic": float(aic) if np.isfinite(aic) else aic,
        "ss_res": float(ss_res), "ss_tot": float(ss_tot),
        "n": int(n)
    }, "ok"

def accept_logistic(lin_r2, lin_aic, log_fit):
    """
    Accept logistic over linear if:
      - R2 improves by LOGISTIC_R2_EPS, or
      - AIC improves (if enabled)
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
# LOGBOOK READING: ROBUST BLOCK SELECTION
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
    xls = pd.ExcelFile(log_xlsx)
    sip_set = set([m for m in sip_measurements_norm if not _is_blank(m)])

    best = None  # (score, sh, blk, id_col, sat_col, mw_col, por_col, vol_col, overlap)

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
    if not results_dir.exists():
        return
    for p in results_dir.iterdir():
        try:
            if p.is_file() or p.is_symlink():
                p.unlink()
        except Exception:
            pass


# =============================================================================
# GNUPLOT UTILITIES
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

def write_plot_xy(results_dir: Path,
                  csv_name: str,
                  out_png: str,
                  source_label: str,
                  title: str,
                  xlabel: str,
                  ylabel: str,
                  x_expr: str,
                  y_expr: str,
                  y_lin_col: str,
                  y_poly_col: str|None,
                  y_log4_col: str|None,
                  legend_lin: str,
                  legend_poly: str|None,
                  legend_log4: str|None):
    """
    Generic single-panel plot for x-y with optional overlays from CSV columns.
    Uses column("name") so headers robust.
    """
    poly_clause = ""
    if y_poly_col:
        poly_clause = f", \\\n     '{csv_name}' using ({x_expr}):({y_poly_col}) with lines dt 1 lw 4 title '{legend_poly}'"
    log4_clause = ""
    if y_log4_col:
        log4_clause = f", \\\n     '{csv_name}' using ({x_expr}):({y_log4_col}) with lines dt 1 lw 4 title '{legend_log4}'"

    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1400,900 enhanced font 'Arial,26'
set output '{out_png}'

set grid
set tics out
set key top right
set title '{source_label}: {title}'
set xlabel '{xlabel}'
set ylabel '{ylabel}'

plot '{csv_name}' using ({x_expr}):({y_expr}) with points pt 7 ps 1.8 title 'data', \\
     '{csv_name}' using ({x_expr}):({y_lin_col}) with lines dt 2 lw 4 title '{legend_lin}'{poly_clause}{log4_clause}
"""
    (results_dir / (Path(out_png).stem + ".gp")).write_text(gp)


# =============================================================================
# FIT RESULT STRUCTURES
# =============================================================================

@dataclass
class FitResult:
    kind: str            # "linear" | "poly" | "logistic4"
    x_name: str
    y_name: str
    normalize: bool
    space: str           # "LOGLOG" or "REAL"
    n_points: int
    slope: float = np.nan
    intercept: float = np.nan
    r2: float = np.nan
    p_value: float = np.nan
    aic: float = np.nan
    poly_deg: int = 0
    poly_r2: float = np.nan
    log4_L: float = np.nan
    log4_k: float = np.nan
    log4_x0: float = np.nan
    log4_c: float = np.nan
    log4_r2: float = np.nan
    log4_p: float = np.nan
    log4_aic: float = np.nan
    log4_status: str = ""
    log4_accepted: bool = False
    log4_reason: str = ""


# =============================================================================
# CORE: RUN ONE MODE (NORM or UNORM)
# =============================================================================

def apply_mode_and_build_columns(df_in: pd.DataFrame,
                                 normalize: bool,
                                 ici_mode: str):
    """
    Build a dataframe with all real values and logs needed for this mode.

    For both modes:
      - x_real = S
      - x_log  = logS
      - y_real candidates:
          rho, sigma_imag
      - y_log candidates:
          logRho, logSig

    For NORM mode:
      - choose reference point as max S
      - compute RI, ICI, and VN target (always ICI_vn = sig_ref/sig)

    For UNORM mode:
      - no RI/ICI indices needed for the "main" unnormalized conductivity relation.
      - still provides logRho, logSig, etc. for UNORM plots.
    """
    df = df_in.copy()

    # Basic real/log columns
    df["logS"] = df["S"].apply(safe_log10)
    df["logRho"] = df["rho"].apply(safe_log10)
    df["logSig"] = df["sigma_imag"].apply(safe_log10)

    # Defaults for index/target columns (may remain NaN for UNORM)
    for c in [
        "rho_ref", "sig_ref", "logRho_ref", "logSig_ref",
        "RI", "ICI_goal", "ICI_vn",
        "logRI_goal", "logICI_goal",
        "logRI_vn", "logICI_vn"
    ]:
        if c not in df.columns:
            df[c] = np.nan

    if normalize:
        # Reference at max S
        iref = int(np.nanargmax(df["S"].to_numpy(float)))
        rho_ref = float(df.iloc[iref]["rho"])
        sig_ref = float(df.iloc[iref]["sigma_imag"])
        logRho_ref = safe_log10(rho_ref)
        logSig_ref = safe_log10(sig_ref)

        df["rho_ref"] = rho_ref
        df["sig_ref"] = sig_ref
        df["logRho_ref"] = logRho_ref
        df["logSig_ref"] = logSig_ref

        # RI always rho/rho_ref
        df["RI"] = df["rho"] / rho_ref
        df["logRI_goal"] = df["logRho"] - logRho_ref  # log(rho/rho_ref)

        # GOAL/INDEX ICI follows --ici-mode
        if ici_mode == "sig_over_sigref":
            df["ICI_goal"] = df["sigma_imag"] / sig_ref
            df["logICI_goal"] = df["logSig"] - logSig_ref
        else:
            df["ICI_goal"] = sig_ref / df["sigma_imag"]
            df["logICI_goal"] = logSig_ref - df["logSig"]

        # TARGET/VN always uses sig_ref/sig (downtrend)
        df["ICI_vn"] = sig_ref / df["sigma_imag"]
        df["logRI_vn"] = df["logRI_goal"]
        df["logICI_vn"] = logSig_ref - df["logSig"]
    else:
        # UNORM: keep rho/sigma relations; indices left NaN
        pass

    return df


def fit_all_methods(x, y, do_poly: bool, poly_deg: int,
                    do_log4: bool, logistic_fast: bool, logistic_max_seconds: float,
                    logistic_min_points: int):
    """
    Return:
      - linear fit dict
      - poly coeffs and r2 (optional)
      - logistic4 accepted fit dict (optional) + status/reason
    """
    # Linear
    slope, intercept, r2, pval, ss_res, aic = linear_fit(x, y)
    lin = {"slope": slope, "intercept": intercept, "r2": r2, "p": pval, "aic": aic}

    # Poly
    poly = {"deg": 0, "r2": np.nan, "coeffs": None}
    if do_poly and poly_deg >= 2:
        coeffs, _, r2p, _ = poly_fit_predict(x, y, poly_deg)
        poly = {"deg": poly_deg, "r2": r2p, "coeffs": coeffs}

    # Logistic4
    log4 = None
    log4_status = "skipped"
    log4_acc = None
    log4_reason = "n/a"

    if do_log4 and len(np.asarray(x)) >= logistic_min_points:
        log4, log4_status = logistic_fit_noscipy(
            x, y, fast=logistic_fast,
            max_seconds=(logistic_max_seconds if logistic_max_seconds and logistic_max_seconds > 0 else 0.0)
        )
        if log4 is not None:
            ok, reason = accept_logistic(lin["r2"], lin["aic"], log4)
            log4_reason = reason
            if ok:
                log4_acc = log4
    elif do_log4:
        log4_status = "not_enough_points"

    return lin, poly, log4_acc, log4_status, log4_reason


def add_fit_columns(df: pd.DataFrame,
                    x_col: str, y_col: str,
                    y_linfit_col: str,
                    y_polyfit_col: str|None,
                    y_log4fit_col: str|None,
                    lin: dict,
                    poly: dict,
                    log4_acc: dict|None):
    """
    Create fitted curve columns (on the original df rows where x/y finite).
    """
    m = np.isfinite(df[x_col]) & np.isfinite(df[y_col])
    xx = df.loc[m, x_col].to_numpy(float)

    # linear
    df.loc[m, y_linfit_col] = lin["slope"]*xx + lin["intercept"]

    # poly
    if y_polyfit_col and poly.get("coeffs", None) is not None:
        df.loc[m, y_polyfit_col] = np.polyval(poly["coeffs"], xx)

    # logistic4
    if y_log4fit_col and log4_acc is not None:
        L = log4_acc["L"]; k = log4_acc["k"]; x0 = log4_acc["x0"]; c = log4_acc["c"]
        df.loc[m, y_log4fit_col] = logistic4(xx, L, k, x0, c)

    return df


# =============================================================================
# MAIN PROCESSING FOR ONE RUN FOLDER
# =============================================================================

def process_run(run_root: Path,
                target_freq: float,
                sheet_filter: str|None,
                dry_run: bool,
                min_points: int,
                max_df: float,
                do_logistic_legacy: bool,
                logistic_min_points: int,
                logistic_fast: bool,
                logistic_max_seconds: float,
                gnuplot_timeout: int,
                clean_results: bool,
                ici_mode: str,
                both_modes: bool,
                poly_deg: int,
                poly_raw: bool,
                do_logistic4: bool):
    """
    Per run:
      - locate files
      - read SIP data
      - map S from logbook
      - filter
      - run BOTH modes if requested
      - write per-mode CSVs, plots, and per-run Excel
      - return rows for global summary (one row per mode)
    """
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return []

    source_label = pretty_source_name(run_root.name)
    results_dir.mkdir(parents=True, exist_ok=True)
    if clean_results:
        clean_results_dir(results_dir)

    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df, sheet_filter=sheet_filter)
    tag = tag_freq(target_freq)

    if sip_df.empty:
        return [{
            "run": str(run_root), "source": source_label, "mode": "N/A",
            "status": "FAIL", "reason": "No SIP sheets readable / target freq missing", "points": 0
        }]

    sat_map, used_sheet, overlap_cnt, block_rows = read_logbook_bestblock_by_overlap(
        log_xlsx, sip_df["measurement_norm"].tolist()
    )

    sip_df["S"] = sip_df["measurement_norm"].map(sat_map).astype(float)

    # Basic filter for real-valued fits:
    # S, rho, sigma must be finite and positive
    m_real = np.isfinite(sip_df["S"]) & np.isfinite(sip_df["rho"]) & np.isfinite(sip_df["sigma_imag"]) & (sip_df["S"] > 0) & (sip_df["rho"] > 0) & (sip_df["sigma_imag"] > 0)
    dropped = sip_df.loc[~m_real, ["measurement", "S", "rho", "sigma_imag"]].copy()
    if len(dropped) > 0:
        dropped.to_csv(results_dir / f"dropped_rows_{tag}Hz.csv", index=False)
    sip_df = sip_df.loc[m_real].copy()

    if len(sip_df) < min_points:
        return [{
            "run": str(run_root), "source": source_label, "mode": "N/A",
            "status": "FAIL",
            "reason": f"Not enough matched points after real filter (need >= {min_points})",
            "points": int(len(sip_df)),
            "logbook_sheet": used_sheet,
            "logbook_overlap": overlap_cnt,
            "logbook_block_rows": block_rows,
        }]

    # Debug dump before any mode logic
    sip_df.to_csv(results_dir / f"debug_before_mode_{tag}Hz.csv", index=False)

    modes = [True] if not both_modes else [True, False]  # True=NORM, False=UNORM
    mode_rows = []

    # Collect per-run excel sheets
    excel_path = results_dir / f"all_values_{tag}Hz.xlsx"
    excel_sheets = {}

    # Collect per-run fit table per mode
    fits_by_mode = {True: [], False: []}

    for normalize in modes:
        df = apply_mode_and_build_columns(sip_df, normalize=normalize, ici_mode=ici_mode)

        # Define per-mode "spaces"
        # REAL space: x=S
        # LOGLOG space: x=logS, and y=log10(y_real)
        # For conductivity goal: sigma_imag vs S is always present.
        # For NORM mode, RI/ICI goal/VN are available (log space only; real RI/ICI also possible).
        # We'll produce both.

        # Prepare output CSV name
        csv_name = f"joined__{mode_suffix(normalize)}_{tag}Hz.csv"
        csv_path = results_dir / csv_name

        # Add columns for fitted curves (we will fill as we go)
        # REAL-space fits
        for c in ["rho_lin_real","rho_poly_real","rho_log4_real",
                  "sig_lin_real","sig_poly_real","sig_log4_real"]:
            df[c] = np.nan

        # LOGLOG fits (on logs)
        for c in ["logRho_lin","logRho_poly","logRho_log4",
                  "logSig_lin","logSig_poly","logSig_log4"]:
            df[c] = np.nan

        # NORM-only index/target fitted curves (log space)
        if normalize:
            for c in ["logRI_goal_lin","logRI_goal_poly","logRI_goal_log4",
                      "logICI_goal_lin","logICI_goal_poly","logICI_goal_log4",
                      "logRI_vn_lin","logRI_vn_poly","logRI_vn_log4",
                      "logICI_vn_lin","logICI_vn_poly","logICI_vn_log4"]:
                df[c] = np.nan

        # ---------- REAL SPACE FITS (rho vs S, sigma vs S) ----------
        xR = df["S"].to_numpy(float)

        # rho vs S (real-real)
        y_rho = df["rho"].to_numpy(float)
        lin, poly, log4_acc, log4_status, log4_reason = fit_all_methods(
            xR, y_rho,
            do_poly=(poly_deg >= 2),
            poly_deg=poly_deg,
            do_log4=do_logistic4,
            logistic_fast=logistic_fast,
            logistic_max_seconds=logistic_max_seconds,
            logistic_min_points=logistic_min_points
        )
        df = add_fit_columns(df, "S", "rho", "rho_lin_real",
                             "rho_poly_real" if (poly_deg >= 2) else None,
                             "rho_log4_real" if do_logistic4 else None,
                             lin, poly, log4_acc)

        fr = FitResult(kind="linear", x_name="S", y_name="rho", normalize=normalize, space="REAL",
                       n_points=int(np.isfinite(xR).sum()), slope=lin["slope"], intercept=lin["intercept"],
                       r2=lin["r2"], p_value=lin["p"], aic=lin["aic"])
        fits_by_mode[normalize].append(asdict(fr))
        if poly_deg >= 2:
            fp = FitResult(kind="poly", x_name="S", y_name="rho", normalize=normalize, space="REAL",
                           n_points=int(np.isfinite(xR).sum()), poly_deg=poly_deg, poly_r2=poly["r2"])
            fits_by_mode[normalize].append(asdict(fp))
        if do_logistic4:
            fl = FitResult(kind="logistic4", x_name="S", y_name="rho", normalize=normalize, space="REAL",
                           n_points=int(np.isfinite(xR).sum()),
                           log4_status=log4_status, log4_accepted=(log4_acc is not None), log4_reason=log4_reason)
            if log4_acc is not None:
                fl.log4_L = log4_acc["L"]; fl.log4_k = log4_acc["k"]; fl.log4_x0 = log4_acc["x0"]; fl.log4_c = log4_acc["c"]
                fl.log4_r2 = log4_acc["r2"]; fl.log4_p = log4_acc["p_model"]; fl.log4_aic = log4_acc["aic"]
            fits_by_mode[normalize].append(asdict(fl))

        # sigma_imag vs S (real-real)  -- KEY GOAL
        y_sig = df["sigma_imag"].to_numpy(float)
        linS, polyS, log4S_acc, log4S_status, log4S_reason = fit_all_methods(
            xR, y_sig,
            do_poly=(poly_deg >= 2),
            poly_deg=poly_deg,
            do_log4=do_logistic4,
            logistic_fast=logistic_fast,
            logistic_max_seconds=logistic_max_seconds,
            logistic_min_points=logistic_min_points
        )
        df = add_fit_columns(df, "S", "sigma_imag", "sig_lin_real",
                             "sig_poly_real" if (poly_deg >= 2) else None,
                             "sig_log4_real" if do_logistic4 else None,
                             linS, polyS, log4S_acc)

        frs = FitResult(kind="linear", x_name="S", y_name="sigma_imag", normalize=normalize, space="REAL",
                        n_points=int(np.isfinite(xR).sum()), slope=linS["slope"], intercept=linS["intercept"],
                        r2=linS["r2"], p_value=linS["p"], aic=linS["aic"])
        fits_by_mode[normalize].append(asdict(frs))
        if poly_deg >= 2:
            fps = FitResult(kind="poly", x_name="S", y_name="sigma_imag", normalize=normalize, space="REAL",
                            n_points=int(np.isfinite(xR).sum()), poly_deg=poly_deg, poly_r2=polyS["r2"])
            fits_by_mode[normalize].append(asdict(fps))
        if do_logistic4:
            fls = FitResult(kind="logistic4", x_name="S", y_name="sigma_imag", normalize=normalize, space="REAL",
                            n_points=int(np.isfinite(xR).sum()),
                            log4_status=log4S_status, log4_accepted=(log4S_acc is not None), log4_reason=log4S_reason)
            if log4S_acc is not None:
                fls.log4_L = log4S_acc["L"]; fls.log4_k = log4S_acc["k"]; fls.log4_x0 = log4S_acc["x0"]; fls.log4_c = log4S_acc["c"]
                fls.log4_r2 = log4S_acc["r2"]; fls.log4_p = log4S_acc["p_model"]; fls.log4_aic = log4S_acc["aic"]
            fits_by_mode[normalize].append(asdict(fls))

        # ---------- LOG-LOG FITS (log10(rho) vs log10(S), log10(sig) vs log10(S)) ----------
        # Need finite logs
        mlog = np.isfinite(df["logS"]) & np.isfinite(df["logRho"]) & np.isfinite(df["logSig"])
        dfl = df.loc[mlog].copy()

        xL = dfl["logS"].to_numpy(float)

        # logRho vs logS
        yLR = dfl["logRho"].to_numpy(float)
        linLR, polyLR, log4LR_acc, log4LR_status, log4LR_reason = fit_all_methods(
            xL, yLR,
            do_poly=(poly_deg >= 2 and poly_raw),  # keep your "poly raw" behavior
            poly_deg=poly_deg,
            do_log4=do_logistic4,
            logistic_fast=logistic_fast,
            logistic_max_seconds=logistic_max_seconds,
            logistic_min_points=logistic_min_points
        )
        # write fits back into df only on mlog
        df.loc[mlog, :] = add_fit_columns(df.loc[mlog, :], "logS", "logRho", "logRho_lin",
                                          "logRho_poly" if (poly_deg >= 2 and poly_raw) else None,
                                          "logRho_log4" if do_logistic4 else None,
                                          linLR, polyLR, log4LR_acc)

        # logSig vs logS
        yLS = dfl["logSig"].to_numpy(float)
        linLS, polyLS, log4LS_acc, log4LS_status, log4LS_reason = fit_all_methods(
            xL, yLS,
            do_poly=(poly_deg >= 2),  # always allow poly here if requested
            poly_deg=poly_deg,
            do_log4=do_logistic4,
            logistic_fast=logistic_fast,
            logistic_max_seconds=logistic_max_seconds,
            logistic_min_points=logistic_min_points
        )
        df.loc[mlog, :] = add_fit_columns(df.loc[mlog, :], "logS", "logSig", "logSig_lin",
                                          "logSig_poly" if (poly_deg >= 2) else None,
                                          "logSig_log4" if do_logistic4 else None,
                                          linLS, polyLS, log4LS_acc)

        # Record key log-log conductivity fits in fit table
        fits_by_mode[normalize].append(asdict(FitResult(
            kind="linear", x_name="logS", y_name="logSig", normalize=normalize, space="LOGLOG",
            n_points=int(len(xL)), slope=linLS["slope"], intercept=linLS["intercept"], r2=linLS["r2"], p_value=linLS["p"], aic=linLS["aic"]
        )))
        if poly_deg >= 2:
            fits_by_mode[normalize].append(asdict(FitResult(
                kind="poly", x_name="logS", y_name="logSig", normalize=normalize, space="LOGLOG",
                n_points=int(len(xL)), poly_deg=poly_deg, poly_r2=polyLS["r2"]
            )))
        if do_logistic4:
            fl = FitResult(kind="logistic4", x_name="logS", y_name="logSig", normalize=normalize, space="LOGLOG",
                           n_points=int(len(xL)), log4_status=log4LS_status, log4_accepted=(log4LS_acc is not None), log4_reason=log4LS_reason)
            if log4LS_acc is not None:
                fl.log4_L = log4LS_acc["L"]; fl.log4_k = log4LS_acc["k"]; fl.log4_x0 = log4LS_acc["x0"]; fl.log4_c = log4LS_acc["c"]
                fl.log4_r2 = log4LS_acc["r2"]; fl.log4_p = log4LS_acc["p_model"]; fl.log4_aic = log4LS_acc["aic"]
            fits_by_mode[normalize].append(asdict(fl))

        # ---------- NORM-only: GOAL + VN target in log space with linear/poly/log4 ----------
        # We keep your original "goal through-origin" concept too, but we also add intercept-based fits here.
        if normalize:
            # GOAL through-origin (original style)
            a_ri, r2_ri, _ = linear_fit_through_origin(dfl["logS"], dfl["logRI_goal"])
            a_ici, r2_ici, _ = linear_fit_through_origin(dfl["logS"], dfl["logICI_goal"])

            # Also do full linear/poly/log4 on those curves (for completeness, log space)
            # logRI_goal vs logS
            y = dfl["logRI_goal"].to_numpy(float)
            linG, polyG, log4G_acc, log4G_status, log4G_reason = fit_all_methods(
                xL, y,
                do_poly=(poly_deg >= 2),
                poly_deg=poly_deg,
                do_log4=do_logistic4,
                logistic_fast=logistic_fast,
                logistic_max_seconds=logistic_max_seconds,
                logistic_min_points=logistic_min_points
            )
            df.loc[mlog, :] = add_fit_columns(df.loc[mlog, :], "logS", "logRI_goal", "logRI_goal_lin",
                                              "logRI_goal_poly" if (poly_deg >= 2) else None,
                                              "logRI_goal_log4" if do_logistic4 else None,
                                              linG, polyG, log4G_acc)

            # logICI_goal vs logS
            y = dfl["logICI_goal"].to_numpy(float)
            linH, polyH, log4H_acc, log4H_status, log4H_reason = fit_all_methods(
                xL, y,
                do_poly=(poly_deg >= 2),
                poly_deg=poly_deg,
                do_log4=do_logistic4,
                logistic_fast=logistic_fast,
                logistic_max_seconds=logistic_max_seconds,
                logistic_min_points=logistic_min_points
            )
            df.loc[mlog, :] = add_fit_columns(df.loc[mlog, :], "logS", "logICI_goal", "logICI_goal_lin",
                                              "logICI_goal_poly" if (poly_deg >= 2) else None,
                                              "logICI_goal_log4" if do_logistic4 else None,
                                              linH, polyH, log4H_acc)

            # VN target (always downward)
            y = dfl["logRI_vn"].to_numpy(float)
            linVN1, polyVN1, log4VN1_acc, _, _ = fit_all_methods(
                xL, y, do_poly=(poly_deg >= 2), poly_deg=poly_deg,
                do_log4=do_logistic4, logistic_fast=logistic_fast,
                logistic_max_seconds=logistic_max_seconds, logistic_min_points=logistic_min_points
            )
            df.loc[mlog, :] = add_fit_columns(df.loc[mlog, :], "logS", "logRI_vn", "logRI_vn_lin",
                                              "logRI_vn_poly" if (poly_deg >= 2) else None,
                                              "logRI_vn_log4" if do_logistic4 else None,
                                              linVN1, polyVN1, log4VN1_acc)

            y = dfl["logICI_vn"].to_numpy(float)
            linVN2, polyVN2, log4VN2_acc, _, _ = fit_all_methods(
                xL, y, do_poly=(poly_deg >= 2), poly_deg=poly_deg,
                do_log4=do_logistic4, logistic_fast=logistic_fast,
                logistic_max_seconds=logistic_max_seconds, logistic_min_points=logistic_min_points
            )
            df.loc[mlog, :] = add_fit_columns(df.loc[mlog, :], "logS", "logICI_vn", "logICI_vn_lin",
                                              "logICI_vn_poly" if (poly_deg >= 2) else None,
                                              "logICI_vn_log4" if do_logistic4 else None,
                                              linVN2, polyVN2, log4VN2_acc)

        # ---------- Write per-mode CSV ----------
        df.to_csv(csv_path, index=False)

        # ---------- Per-mode plots ----------
        # We always create sigma vs S plots (REAL + LOGLOG) for BOTH modes (your final goal).
        # REAL sigma vs S
        out_real_sig = fname("cond_vs_S", normalize, False, tag, extra=("LINEAR" + (f"+POLYdeg{poly_deg}" if poly_deg>=2 else "") + ("+LOGISTIC4" if do_logistic4 else "")))
        write_plot_xy(
            results_dir, csv_path.name, out_real_sig, source_label,
            title=f"Conductivity vs Saturation (REAL) [{mode_suffix(normalize)}]",
            xlabel="Saturation (S)", ylabel="Imag Conductivity (sigma_imag)",
            x_expr='column("S")', y_expr='column("sigma_imag")',
            y_lin_col='column("sig_lin_real")',
            y_poly_col=('column("sig_poly_real")' if poly_deg>=2 else None),
            y_log4_col=('column("sig_log4_real")' if do_logistic4 else None),
            legend_lin="linear",
            legend_poly=(f"poly deg {poly_deg}" if poly_deg>=2 else None),
            legend_log4=("logistic4 (accepted)" if do_logistic4 else None)
        )

        # LOGLOG sigma vs S
        out_log_sig = fname("cond_vs_S", normalize, True, tag, extra=("LINEAR" + (f"+POLYdeg{poly_deg}" if poly_deg>=2 else "") + ("+LOGISTIC4" if do_logistic4 else "")))
        write_plot_xy(
            results_dir, csv_path.name, out_log_sig, source_label,
            title=f"log10(sigma_imag) vs log10(S) [{mode_suffix(normalize)}]",
            xlabel="log10(S)", ylabel="log10(sigma_imag)",
            x_expr='column("logS")', y_expr='column("logSig")',
            y_lin_col='column("logSig_lin")',
            y_poly_col=('column("logSig_poly")' if poly_deg>=2 else None),
            y_log4_col=('column("logSig_log4")' if do_logistic4 else None),
            legend_lin="linear",
            legend_poly=(f"poly deg {poly_deg}" if poly_deg>=2 else None),
            legend_log4=("logistic4 (accepted)" if do_logistic4 else None)
        )

        # Also rho vs S plots (REAL + LOGLOG) for BOTH modes
        out_real_rho = fname("rho_vs_S", normalize, False, tag, extra=("LINEAR" + (f"+POLYdeg{poly_deg}" if poly_deg>=2 else "") + ("+LOGISTIC4" if do_logistic4 else "")))
        write_plot_xy(
            results_dir, csv_path.name, out_real_rho, source_label,
            title=f"Resistivity vs Saturation (REAL) [{mode_suffix(normalize)}]",
            xlabel="Saturation (S)", ylabel="Resistivity (rho)",
            x_expr='column("S")', y_expr='column("rho")',
            y_lin_col='column("rho_lin_real")',
            y_poly_col=('column("rho_poly_real")' if poly_deg>=2 else None),
            y_log4_col=('column("rho_log4_real")' if do_logistic4 else None),
            legend_lin="linear",
            legend_poly=(f"poly deg {poly_deg}" if poly_deg>=2 else None),
            legend_log4=("logistic4 (accepted)" if do_logistic4 else None)
        )

        out_log_rho = fname("rho_vs_S", normalize, True, tag, extra=("LINEAR" + (f"+POLYdeg{poly_deg}" if (poly_deg>=2 and poly_raw) else "") + ("+LOGISTIC4" if do_logistic4 else "")))
        write_plot_xy(
            results_dir, csv_path.name, out_log_rho, source_label,
            title=f"log10(rho) vs log10(S) [{mode_suffix(normalize)}]",
            xlabel="log10(S)", ylabel="log10(rho)",
            x_expr='column("logS")', y_expr='column("logRho")',
            y_lin_col='column("logRho_lin")',
            y_poly_col=('column("logRho_poly")' if (poly_deg>=2 and poly_raw) else None),
            y_log4_col=('column("logRho_log4")' if do_logistic4 else None),
            legend_lin="linear",
            legend_poly=(f"poly deg {poly_deg}" if (poly_deg>=2 and poly_raw) else None),
            legend_log4=("logistic4 (accepted)" if do_logistic4 else None)
        )

        # NORM-only: VN target and GOAL index plots in LOGLOG space
        if normalize:
            # VN target: logRI_vn vs logS and logICI_vn vs logS (separate files)
            out_vn_ri = fname("targetVN_RI", True, True, tag, extra=("LINEAR" + (f"+POLYdeg{poly_deg}" if poly_deg>=2 else "") + ("+LOGISTIC4" if do_logistic4 else "")))
            write_plot_xy(
                results_dir, csv_path.name, out_vn_ri, source_label,
                title="TARGET/VN: logRI_vn vs logS",
                xlabel="log10(S)", ylabel="log10(Resistivity Index)",
                x_expr='column("logS")', y_expr='column("logRI_vn")',
                y_lin_col='column("logRI_vn_lin")',
                y_poly_col=('column("logRI_vn_poly")' if poly_deg>=2 else None),
                y_log4_col=('column("logRI_vn_log4")' if do_logistic4 else None),
                legend_lin="linear",
                legend_poly=(f"poly deg {poly_deg}" if poly_deg>=2 else None),
                legend_log4=("logistic4 (accepted)" if do_logistic4 else None)
            )

            out_vn_ici = fname("targetVN_ICI", True, True, tag, extra=("LINEAR" + (f"+POLYdeg{poly_deg}" if poly_deg>=2 else "") + ("+LOGISTIC4" if do_logistic4 else "")))
            write_plot_xy(
                results_dir, csv_path.name, out_vn_ici, source_label,
                title="TARGET/VN: logICI_vn vs logS",
                xlabel="log10(S)", ylabel="log10(Imag Conductivity Index)",
                x_expr='column("logS")', y_expr='column("logICI_vn")',
                y_lin_col='column("logICI_vn_lin")',
                y_poly_col=('column("logICI_vn_poly")' if poly_deg>=2 else None),
                y_log4_col=('column("logICI_vn_log4")' if do_logistic4 else None),
                legend_lin="linear",
                legend_poly=(f"poly deg {poly_deg}" if poly_deg>=2 else None),
                legend_log4=("logistic4 (accepted)" if do_logistic4 else None)
            )

        # ---------- Run gnuplot unless dry-run ----------
        if not dry_run:
            # run the .gp files we generated (they sit next to the png with same stem)
            for png in [out_real_sig, out_log_sig, out_real_rho, out_log_rho]:
                gp = Path(png).stem + ".gp"
                _run_gnuplot(gp, cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)

            if normalize:
                for png in [out_vn_ri, out_vn_ici]:
                    gp = Path(png).stem + ".gp"
                    _run_gnuplot(gp, cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)

        # ---------- Per-run Excel sheets ----------
        # Split into "LOGLOG" and "REAL" data views for this mode.
        # (We keep them as separate sheets to make browsing easy.)
        cols_base = ["measurement","measurement_norm","freq_hz","S","rho","sigma_imag","logS","logRho","logSig",
                     "rho_ref","sig_ref","logRho_ref","logSig_ref","RI","ICI_goal","ICI_vn","logRI_goal","logICI_goal","logRI_vn","logICI_vn"]
        cols_base = [c for c in cols_base if c in df.columns]

        cols_real_fits = ["rho_lin_real","rho_poly_real","rho_log4_real","sig_lin_real","sig_poly_real","sig_log4_real"]
        cols_log_fits  = ["logRho_lin","logRho_poly","logRho_log4","logSig_lin","logSig_poly","logSig_log4"]
        cols_norm_more = []
        if normalize:
            cols_norm_more = ["logRI_vn_lin","logRI_vn_poly","logRI_vn_log4","logICI_vn_lin","logICI_vn_poly","logICI_vn_log4",
                              "logRI_goal_lin","logRI_goal_poly","logRI_goal_log4","logICI_goal_lin","logICI_goal_poly","logICI_goal_log4"]

        # REAL sheet: show raw + real-fit columns
        sheet_real = f"data__{mode_suffix(normalize)}__REAL"
        excel_sheets[sheet_real] = df[cols_base + cols_real_fits + cols_norm_more].copy()

        # LOGLOG sheet: show logs + log-fit columns
        sheet_log = f"data__{mode_suffix(normalize)}__LOGLOG"
        excel_sheets[sheet_log] = df[cols_base + cols_log_fits + cols_norm_more].copy()

        # ---------- Per-mode summary row ----------
        mode_rows.append({
            "run": str(run_root),
            "source": source_label,
            "mode": mode_suffix(normalize),
            "status": "OK",
            "points": int(len(df)),
            "logbook_sheet": used_sheet,
            "logbook_overlap": overlap_cnt,
            "logbook_block_rows": block_rows,
            "ici_mode": ici_mode,
            "tagHz": tag,
            "joined_csv": str(csv_path),
            "per_run_excel": str(excel_path),
        })

    # Write per-run Excel (once, after both modes done)
    # Also write fit summaries as separate sheets.
    fits_norm_df = pd.DataFrame(fits_by_mode[True]) if fits_by_mode.get(True) else pd.DataFrame()
    fits_un_df   = pd.DataFrame(fits_by_mode[False]) if fits_by_mode.get(False) else pd.DataFrame()

    with pd.ExcelWriter(excel_path, engine="openpyxl") as xw:
        # Data sheets
        for sh, dfx in excel_sheets.items():
            # Excel sheet names max length 31
            sh2 = sh[:31]
            dfx.to_excel(xw, sheet_name=sh2, index=False)

        # Fit summary sheets
        if not fits_norm_df.empty:
            fits_norm_df.to_excel(xw, sheet_name="fits__NORM", index=False)
        if not fits_un_df.empty:
            fits_un_df.to_excel(xw, sheet_name="fits__UNORM", index=False)

    return mode_rows


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
# MAIN CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("top", help="Top folder containing run subfolders")
    ap.add_argument("--freq", type=float, default=0.01, help="Target frequency in Hz")
    ap.add_argument("--sample", type=str, default=None, help="Optional substring filter on SIP sheet names")
    ap.add_argument("--dry-run", action="store_true", help="Write scripts/CSVs but do not run gnuplot")
    ap.add_argument("--min-points", type=int, default=4, help="Minimum matched points to accept a run")
    ap.add_argument("--max-df", type=float, default=0.0, help="Max allowed |f_found - f_target| (0 disables)")

    ap.add_argument("--logistic-min-points", type=int, default=6, help="Min points needed to attempt logistic fits")
    ap.add_argument("--logistic-fast", action="store_true", help="Use SciPy-free fast logistic search")
    ap.add_argument("--logistic-max-seconds", type=float, default=8.0, help="Time cap for logistic fit per run")
    ap.add_argument("--gnuplot-timeout", type=int, default=30, help="Timeout seconds per gnuplot call")

    ap.add_argument("--clean-results", action="store_true",
                    help="Delete files inside each run's Results/ before writing new outputs.")

    ap.add_argument("--ici-mode", type=str, default="sigref_over_sig",
                    choices=["sigref_over_sig","sig_over_sigref"],
                    help=("NORM mode GOAL/INDEX ICI definition. VN target is fixed to sigref/sig "
                          "so that target trends downward."))

    ap.add_argument("--poly-deg", type=int, default=0,
                    help="Polynomial degree for extra fits/plots (>=2). 0 disables.")
    ap.add_argument("--poly-raw", action="store_true",
                    help="Apply polynomial overlay to RAW logRho vs logS (LOGLOG) (requires --poly-deg>=2).")

    # NEW: run both normalization modes in one invocation
    ap.add_argument("--both-modes", action="store_true",
                    help="Run both NORM and UNORM modes (all methods/spaces) in one pass per folder.")

    # NEW: enable logistic4 broadly (REAL and LOGLOG across modes)
    ap.add_argument("--logistic4", action="store_true",
                    help="Enable 4-parameter logistic fits broadly in REAL and LOGLOG spaces for rho and sigma.")

    args = ap.parse_args()

    base = Path(args.top).expanduser().resolve()
    if not base.exists():
        print(f"Not found: {base}")
        raise SystemExit(2)

    roots = discover_run_roots_fast(base)
    print(f"Found {len(roots)} run folder(s).", flush=True)

    rows_all = []
    t0 = time.time()

    for i, run in enumerate(roots, 1):
        print(f"[{i}/{len(roots)}] START {run}", flush=True)
        ti = time.time()
        try:
            mode_rows = process_run(
                run_root=run,
                target_freq=args.freq,
                sheet_filter=args.sample,
                dry_run=args.dry_run,
                min_points=args.min_points,
                max_df=args.max_df,
                do_logistic_legacy=True,               # kept placeholder (legacy v27 logistic-only plots are superseded by broad logistic4)
                logistic_min_points=args.logistic_min_points,
                logistic_fast=args.logistic_fast,
                logistic_max_seconds=args.logistic_max_seconds,
                gnuplot_timeout=args.gnuplot_timeout,
                clean_results=args.clean_results,
                ici_mode=args.ici_mode,
                both_modes=args.both_modes,
                poly_deg=args.poly_deg,
                poly_raw=args.poly_raw,
                do_logistic4=args.logistic4
            )
            rows_all.extend(mode_rows)
            dt = time.time() - ti
            # compact print
            if mode_rows:
                ok_modes = ", ".join([r.get("mode","?") for r in mode_rows if r.get("status") == "OK"])
                pts = mode_rows[0].get("points", 0)
                print(f"[{i}/{len(roots)}] OK modes=[{ok_modes}] points={pts} ({dt:.1f}s)", flush=True)
            else:
                print(f"[{i}/{len(roots)}] SKIP (missing inputs) ({dt:.1f}s)", flush=True)

        except subprocess.TimeoutExpired:
            rows_all.append({"run": str(run), "status": "FAIL", "reason": "gnuplot timeout", "points": 0})
            print(f"[{i}/{len(roots)}] FAIL gnuplot timeout", flush=True)
        except Exception as e:
            rows_all.append({"run": str(run), "status": "FAIL", "reason": str(e), "points": 0})
            print(f"[{i}/{len(roots)}] FAIL {e}", flush=True)

    tag = tag_freq(args.freq)
    summary_csv = base / f"Results_Summary_{tag}Hz.csv"
    summary_xlsx = base / f"Results_Summary_{tag}Hz.xlsx"

    if rows_all:
        df_sum = pd.DataFrame(rows_all)
        df_sum.to_csv(summary_csv, index=False)
        with pd.ExcelWriter(summary_xlsx, engine="openpyxl") as xw:
            df_sum.to_excel(xw, sheet_name="summary", index=False)
        print(f"\nSummary written: {summary_csv}", flush=True)
        print(f"Summary written: {summary_xlsx}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)

if __name__ == "__main__":
    main()
