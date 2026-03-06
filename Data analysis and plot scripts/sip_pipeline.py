#!/usr/bin/env python3
"""
sip_pipeline.py

legacy-style + requested updates:
----------------------------------------------------
1) LOGISTIC4 runs with >=4 points.
2) Restores legacy plots (RAW/INDEX/TARGET/diagnostics/logistic-only):
   - same naming
   - same "feel" (fonts/points/linetypes)
   - RAW/INDEX/TARGET are tall again (1200x1600) so TARGET is not "wide".
3) Puts fit parameters back ON the legacy plots (labels):
   - RAW: slope/intercept/R2/p
   - INDEX: slope/intercept/R2/p
   - TARGET (VN-style through-origin): exponent/R2/SSE
   - Diagnostics: linear + logistic params/status
4) --clean-results no longer deletes old plots:
   - it backups Results/ to Results__backup_<timestamp>/ instead of deleting.
5) Keeps extra plots:
   - REAL and LOGLOG bundles for rho vs S and conductivity vs S
   - BOTH modes: NORM + UNORM (via --both-modes)
   - Poly fits (deg>=2) on the bundle plots and stored in CSV
   - LOG4 fits (real-space and log-space) stored in CSV
6) Adds fit equations for all fits:
   - equations sheet in Excel outputs always
   - optional equation labels in the right margin of plots via --plot-equations
     (default: OFF)

FIXES:
------------------------
A) No more gnuplot spam: write CSVs with na_rep="" (blank instead of 'nan')
B) Legacy logistic plots ALWAYS show logistic4 if it was computed
   (even if "rejected" by AIC/R2 rule). We label accepted/rejected
   instead of hiding the curve.
C) Diagnostics/logistic-only titles/labels show acceptance status + reason.
D) Added missing LOGLOG conductivity logistic fit:
   - logSig_log4 is computed from logSig vs logS
E) Added back old target-result plot as an extra figure:
   - target_raw_exponents_<tag>Hz.png
   - uses linear-with-intercept results (blue/red pair logic)
F) Fixed titles such as VN S-12 to stay plain text (no subscript conversion).

Directory assumptions:
----------------------------------
<run_root>/
  Analysis/         -> SIP output Excel (filename contains "output" usually)
  SIP_Raw_Data/     -> logbook Excel (filename contains "Log Book")
  Results/          -> created and filled with csv + plots + reports + excel
"""

import math
import argparse
import subprocess
import time
import re
import shutil
from datetime import datetime
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
    return s

def tag_from_freq(freq: float) -> str:
    return str(freq).replace(".", "p")

def fmt(x, nd=4):
    if x is None:
        return "nan"
    try:
        x = float(x)
        if not np.isfinite(x):
            return "nan"
        return f"{x:.{nd}g}"
    except Exception:
        return "nan"

def gpeq(s: str) -> str:
    return str(s).replace("\\", "\\\\").replace("'", "\\'")

def eq_linear_str(a, b, x="x", y="y"):
    if not (np.isfinite(a) and np.isfinite(b)):
        return f"{y} = nan"
    sign = "+" if float(b) >= 0 else "-"
    return f"{y} = {fmt(a,6)}*{x} {sign} {fmt(abs(float(b)),6)}"

def eq_origin_str(a, x="x", y="y"):
    if not np.isfinite(a):
        return f"{y} = nan"
    return f"{y} = {fmt(a,6)}*{x}"

def eq_poly_str(coeffs, x="x", y="y"):
    if coeffs is None:
        return f"{y} = nan"
    coeffs = list(coeffs)
    deg = len(coeffs) - 1
    parts = []
    for i, c in enumerate(coeffs):
        p = deg - i
        if not np.isfinite(c):
            continue
        cabs = fmt(abs(float(c)), 6)
        if p == 0:
            term = f"{cabs}"
        elif p == 1:
            term = f"{cabs}*{x}"
        else:
            term = f"{cabs}*{x}^{p}"

        if not parts:
            if float(c) < 0:
                parts.append(f"-{term}")
            else:
                parts.append(term)
        else:
            parts.append((" - " if float(c) < 0 else " + ") + term)

    if not parts:
        return f"{y} = nan"
    return f"{y} = " + "".join(parts)

def eq_log4_str(L, k, x0, c, x="x", y="y"):
    vals = [L, k, x0, c]
    if not all(np.isfinite(v) for v in vals):
        return f"{y} = nan"
    sign_x0 = "-" if float(x0) >= 0 else "+"
    return (
        f"{y} = {fmt(c,6)} + {fmt(L,6)}/"
        f"(1 + exp(-{fmt(k,6)}*({x} {sign_x0} {fmt(abs(float(x0)),6)})))"
    )

def add_equation_row(eq_rows: list, run_root: Path, source_label: str, mode: str,
                     freq: float, plot_name: str, fit_name: str, equation: str,
                     r2=np.nan, p=np.nan, aic=np.nan, accepted="", reason=""):
    eq_rows.append({
        "run_root": str(run_root),
        "source": source_label,
        "mode": mode,
        "freq_target": float(freq),
        "plot": plot_name,
        "fit_name": fit_name,
        "equation": equation,
        "r2": r2,
        "p_value": p,
        "aic": aic,
        "accepted": accepted,
        "reason": reason,
    })

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
# FIT ROUTINES: LINEAR / THROUGH ORIGIN / POLY
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
# LOGISTIC FIT (NO SCIPY) - 4 PARAM
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

def read_sip_output(sip_xlsx: Path, target_freq: float, max_df: float, sheet_filter: str | None):
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
            Mw  = to_float(r.get(mw_col, np.nan))
            phi = to_float(r.get(por_col, np.nan))
            V   = to_float(r.get(vol_col, np.nan))
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

def backup_results_dir(results_dir: Path):
    if not results_dir.exists():
        return
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = results_dir.parent / f"{results_dir.name}__backup_{ts}"
    try:
        if backup.exists():
            shutil.rmtree(backup, ignore_errors=True)
        shutil.move(str(results_dir), str(backup))
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
# LEGACY GNUPLOT WRITERS
# =============================================================================

def write_plot_raw_exponents(results_dir: Path, csv_name: str, out_png: str, freq: float,
                             source_label: str,
                             rho_slope, rho_int, rho_r2, rho_p, rho_eq,
                             sig_slope, sig_int, sig_r2, sig_p, sig_eq,
                             show_eq: bool):
    tag = "{:g}".format(freq)
    eq_top = f"set label 11 '{gpeq(rho_eq)}' at graph 1.02,0.84 left front\n" if show_eq else ""
    eq_bot = f"set label 12 '{gpeq(sig_eq)}' at graph 1.02,0.84 left front\n" if show_eq else ""
    rm = "set rmargin 34" if show_eq else "set rmargin 10"
    un_top = "unset label 11" if show_eq else ""
    un_bot = "unset label 12" if show_eq else ""
    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 2,1

set grid
set key top right
set tics out
{rm}

# --- TOP: logRho vs logS
set title '{gpeq(source_label)}: RAW exponents (log-log) at {tag} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
set label 1 sprintf('linear: a={fmt(rho_slope)}  b={fmt(rho_int)}  R^2={fmt(rho_r2)}  p={fmt(rho_p)}') at graph 0.02,0.94 front
{eq_top}plot '{csv_name}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logRho_lin")) with lines dt 2 lw 3 title 'linear'
unset label 1
{un_top}

# --- BOTTOM: logSig vs logS
set title '{gpeq(source_label)}: RAW exponents (log-log) at {tag} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Conductivity)'
set label 2 sprintf('linear: a={fmt(sig_slope)}  b={fmt(sig_int)}  R^2={fmt(sig_r2)}  p={fmt(sig_p)}') at graph 0.02,0.94 front
{eq_bot}plot '{csv_name}' using (column("logS")):(column("logSig")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logSig_lin")) with lines dt 2 lw 3 title 'linear'
unset label 2
{un_bot}

unset multiplot
"""
    (results_dir / "plot_raw_exponents.gp").write_text(gp)

def write_plot_index_exponents(results_dir: Path, csv_name: str, out_png: str, freq: float,
                               source_label: str,
                               ri_slope, ri_int, ri_r2, ri_p, ri_eq,
                               ici_slope, ici_int, ici_r2, ici_p, ici_eq,
                               show_eq: bool):
    tag = "{:g}".format(freq)
    eq_top = f"set label 11 '{gpeq(ri_eq)}' at graph 1.02,0.84 left front\n" if show_eq else ""
    eq_bot = f"set label 12 '{gpeq(ici_eq)}' at graph 1.02,0.84 left front\n" if show_eq else ""
    rm = "set rmargin 34" if show_eq else "set rmargin 10"
    un_top = "unset label 11" if show_eq else ""
    un_bot = "unset label 12" if show_eq else ""
    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 2,1

set grid
set key top right
set tics out
{rm}

# --- TOP: logRhoIdx vs logS
set title '{gpeq(source_label)}: INDEX exponents at {tag} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(RI) = log10(rho/rho_ref)'
set label 1 sprintf('linear: a={fmt(ri_slope)}  b={fmt(ri_int)}  R^2={fmt(ri_r2)}  p={fmt(ri_p)}') at graph 0.02,0.94 front
{eq_top}plot '{csv_name}' using (column("logS")):(column("logRhoIdx")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logRhoIdx_lin")) with lines dt 2 lw 3 title 'linear'
unset label 1
{un_top}

# --- BOTTOM: logICI (VN-style index) vs logS
set title '{gpeq(source_label)}: INDEX exponents at {tag} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(ICI) index'
set label 2 sprintf('linear: a={fmt(ici_slope)}  b={fmt(ici_int)}  R^2={fmt(ici_r2)}  p={fmt(ici_p)}') at graph 0.02,0.94 front
{eq_bot}plot '{csv_name}' using (column("logS")):(column("logICI_vn")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logICI_lin")) with lines dt 2 lw 3 title 'linear'
unset label 2
{un_bot}

unset multiplot
"""
    (results_dir / "plot_index_exponents.gp").write_text(gp)

def write_plot_target_exponents(results_dir: Path, csv_name: str, out_png: str, freq: float,
                                source_label: str,
                                goal_ri_a, goal_ri_r2, goal_ri_sse, goal_ri_eq,
                                goal_ici_a, goal_ici_r2, goal_ici_sse, goal_ici_eq,
                                show_eq: bool):
    tag = "{:g}".format(freq)
    eq_top = f"set label 11 '{gpeq(goal_ri_eq)}' at graph 1.02,0.84 left front\n" if show_eq else ""
    eq_bot = f"set label 12 '{gpeq(goal_ici_eq)}' at graph 1.02,0.84 left front\n" if show_eq else ""
    rm = "set rmargin 34" if show_eq else "set rmargin 10"
    un_top = "unset label 11" if show_eq else ""
    un_bot = "unset label 12" if show_eq else ""
    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 2,1

set grid
set key top right
set tics out
{rm}

# --- TOP: TARGET RI (VN)
set title '{gpeq(source_label)}: TARGET (VN-style) at {tag} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'logRI (VN) = log10(rho/rho_ref)'
set label 1 sprintf('through-origin: a={fmt(goal_ri_a)}  R^2={fmt(goal_ri_r2)}  SSE={fmt(goal_ri_sse)}') at graph 0.02,0.94 front
{eq_top}plot '{csv_name}' using (column("logS")):(column("logRI_vn")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logRI_goal")) with lines dt 2 lw 3 title 'through-origin'
unset label 1
{un_top}

# --- BOTTOM: TARGET ICI (VN)
set title '{gpeq(source_label)}: TARGET (VN-style) at {tag} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'logICI (VN) = log10(sigref/sig)'
set label 2 sprintf('through-origin: a={fmt(goal_ici_a)}  R^2={fmt(goal_ici_r2)}  SSE={fmt(goal_ici_sse)}') at graph 0.02,0.94 front
{eq_bot}plot '{csv_name}' using (column("logS")):(column("logICI_vn")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logICI_goal")) with lines dt 2 lw 3 title 'through-origin'
unset label 2
{un_bot}

unset multiplot
"""
    (results_dir / "plot_target_exponents.gp").write_text(gp)

def write_plot_target_raw_exponents(results_dir: Path, csv_name: str, out_png: str, freq: float,
                                    source_label: str,
                                    n_exp, n_r2, n_eq,
                                    p_exp, p_r2, p_eq,
                                    show_eq: bool):
    tag = "{:g}".format(freq)
    eq_top = f"set label 11 '{gpeq(n_eq)}' at graph 1.02,0.84 left front\n" if show_eq else ""
    eq_bot = f"set label 12 '{gpeq(p_eq)}' at graph 1.02,0.84 left front\n" if show_eq else ""
    rm = "set rmargin 34" if show_eq else "set rmargin 10"
    un_top = "unset label 11" if show_eq else ""
    un_bot = "unset label 12" if show_eq else ""
    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 2,1

set grid
set key top right
set tics out
{rm}

# --- TOP
set title '{gpeq(source_label)}: Saturation Exponent (Resistivity) at {tag} Hz'
set xlabel 'log(Saturation)'
set ylabel 'log(Resistivity Index)'
set label 1 sprintf('{tag}Hz: n = {fmt(n_exp)}  R^2 = {fmt(n_r2)}') at graph 0.02,0.94 front
{eq_top}plot '{csv_name}' using (column("logS")):(column("logRhoIdx")) with points pt 7 ps 1.8 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logRhoIdx_lin")) with lines dt 2 lw 3 title sprintf('{tag}Hz: n = {fmt(n_exp)}, R^2 = {fmt(n_r2)}')
unset label 1
{un_top}

# --- BOTTOM
set title '{gpeq(source_label)}: Saturation Exponent (Imaginary Conductivity) at {tag} Hz'
set xlabel 'log(Saturation)'
set ylabel 'log(Imaginary Conductivity Index)'
set label 2 sprintf('{tag}Hz: p = {fmt(p_exp)}  R^2 = {fmt(p_r2)}') at graph 0.02,0.94 front
{eq_bot}plot '{csv_name}' using (column("logS")):(column("logICI_vn")) with points pt 7 ps 1.8 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logICI_lin")) with lines dt 2 lw 3 title sprintf('{tag}Hz: p = {fmt(p_exp)}, R^2 = {fmt(p_r2)}')
unset label 2
{un_bot}

unset multiplot
"""
    (results_dir / "plot_target_raw_exponents.gp").write_text(gp)

def write_plot_logistic_diag(results_dir: Path, csv_name: str, out_png: str,
                             freq: float, source_label: str,
                             lin_slope, lin_int, lin_r2, lin_eq,
                             has_log_points: bool,
                             log4_params: dict | None,
                             log4_eq: str,
                             log4_accepted: bool,
                             log4_reason: str,
                             show_eq: bool):
    tag = "{:g}".format(freq)

    extra_curve = ""
    log_label = ""
    log_eq_label = ""

    if has_log_points:
        status_txt = "accepted" if log4_accepted else "rejected"
        extra_curve = (
            ", \\\n     '{csv}' using (column(\"logS\")):(column(\"logRho_logfit\")) "
            "with lines dt 1 lw 4 title 'logistic4 ({status})'"
        ).format(csv=csv_name, status=status_txt)

        if log4_params:
            log_label = (
                f"set label 9 sprintf('log4 ({status_txt}): L={fmt(log4_params.get('L'))}  k={fmt(log4_params.get('k'))}  "
                f"x0={fmt(log4_params.get('x0'))}  c={fmt(log4_params.get('c'))}  "
                f"R^2={fmt(log4_params.get('r2'))}  AIC={fmt(log4_params.get('aic'))}  ({gpeq(log4_reason)})') at graph 0.02,0.84 front\n"
            )
            if show_eq:
                log_eq_label = f"set label 19 '{gpeq(log4_eq)}' at graph 1.02,0.72 left front\n"

    if has_log_points:
        panel3 = f"plot '{csv_name}' using (column(\"logS\")):(column(\"resid_log_rho\")) with points pt 7 ps 1.4 title 'residuals'\n"
        extra3 = ""
    else:
        extra3 = (
            "unset key\n"
            "set label 8 'No logistic fit' at graph 0.5,0.5 center front\n"
            "set xrange [-2:0]\n"
            "set yrange [-1:1]\n"
        )
        panel3 = (
            "plot '-' using 1:2 with points pt 7 ps 0 notitle\n"
            "-2 0\n"
            "0  0\n"
            "e\n"
        )

    rm = "set rmargin 38" if show_eq else "set rmargin 10"
    lin_eq_label = f"set label 11 '{gpeq(lin_eq)}' at graph 1.02,0.90 left front\n" if show_eq else ""
    un11 = "unset label 11" if show_eq else ""
    un19 = "unset label 19" if show_eq else ""

    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1800 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 3,1

set title '{gpeq(source_label)}: Diagnostics: logRho vs logS at {tag} Hz'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
{rm}
set label 1 sprintf('linear: a={fmt(lin_slope)}  b={fmt(lin_int)}  R^2={fmt(lin_r2)}') at graph 0.02,0.92 front
{lin_eq_label}{log_label}{log_eq_label}plot '{csv_name}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logRho_linfit")) with lines dt 2 lw 3 title 'linear'{extra_curve}
unset label 1
{un11}
unset label 9
{un19}

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
{extra3}{panel3}
unset label 8

unset multiplot
"""
    (results_dir / "plot_logistic_diagnostics.gp").write_text(gp)

def write_plot_logistic_only(results_dir: Path, csv_name: str, out_png: str, freq: float,
                             source_label: str, has_log_points: bool,
                             log4_eq: str,
                             log4_accepted: bool, log4_reason: str,
                             show_eq: bool):
    tag = "{:g}".format(freq)

    if has_log_points:
        status_txt = "accepted" if log4_accepted else "rejected"
        rm = "set rmargin 34" if show_eq else "set rmargin 10"
        eq_label = f"set label 1 '{gpeq(log4_eq)}' at graph 1.02,0.92 left front\n" if show_eq else ""
        un1 = "unset label 1" if show_eq else ""
        gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,24'
set output '{out_png}'
set title '{gpeq(source_label)}: Logistic-only (logRho vs logS) at {tag} Hz  [log4: {status_txt} | {gpeq(log4_reason)}]'

set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
{rm}
{eq_label}plot '{csv_name}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using (column("logS")):(column("logRho_logfit")) with lines dt 1 lw 4 title 'logistic4 ({status_txt})'
{un1}
"""
    else:
        gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,24'
set output '{out_png}'
set title '{gpeq(source_label)}: Logistic-only (logRho vs logS) at {tag} Hz'

unset key
set grid
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
set label 1 'No logistic fit' at graph 0.5,0.5 center front
set xrange [-2:0]
set yrange [0:1]
plot '-' using 1:2 with points pt 7 ps 0 notitle
-2 0.5
0  0.5
e
"""
    (results_dir / "plot_logistic_only.gp").write_text(gp)

# =============================================================================
# BUNDLE PLOTS
# =============================================================================

def write_plot_xy_bundle(results_dir: Path, csv_name: str, out_png: str, out_gp: str,
                         source_label: str, freq: float,
                         x_expr: str, y_expr: str,
                         xlab: str, ylab: str,
                         title: str,
                         do_linear: bool, do_poly: bool, poly_deg: int,
                         do_log4: bool,
                         y_lin_col: str, y_poly_col: str, y_log4_col: str,
                         linear_eq: str = "",
                         poly_eq: str = "",
                         log4_eq: str = "",
                         show_eq: bool = False):
    label = "{:g}".format(freq)

    df = None
    try:
        df = pd.read_csv(results_dir / csv_name)
    except Exception:
        df = None

    def col_has_finite(colname: str) -> bool:
        if df is None:
            return True
        if colname not in df.columns:
            return False
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

    label_lines = []
    if show_eq:
        y0 = 0.94
        if do_linear and linear_eq:
            label_lines.append(f"set label 101 '{gpeq(linear_eq)}' at graph 1.02,{y0:.2f} left front")
            y0 -= 0.08
        if do_poly and poly_deg >= 2 and poly_eq:
            label_lines.append(f"set label 102 '{gpeq(poly_eq)}' at graph 1.02,{y0:.2f} left front")
            y0 -= 0.10
        if do_log4 and log4_eq:
            label_lines.append(f"set label 103 '{gpeq(log4_eq)}' at graph 1.02,{y0:.2f} left front")
            y0 -= 0.08

    unset_lines = []
    if show_eq:
        for lid in [101, 102, 103]:
            unset_lines.append(f"unset label {lid}")

    rm = "set rmargin 42" if show_eq else "set rmargin 10"

    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,24'
set output '{out_png}'
set grid
set key top right
set tics out
{rm}

set title '{gpeq(title)} ({label} Hz)'
set xlabel '{gpeq(xlab)}'
set ylabel '{gpeq(ylab)}'

{chr(10).join(label_lines)}
{plot_cmd}
{chr(10).join(unset_lines)}
"""
    (results_dir / out_gp).write_text(gp)

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
# MAIN PROCESSING FOR ONE RUN (ONE MODE)
# =============================================================================

def process_one_mode(run_root: Path, target_freq: float, sheet_filter: str | None,
                     min_points: int, max_df: float,
                     do_log4: bool, logistic_fast: bool, logistic_max_seconds: float,
                     gnuplot_timeout: int,
                     clean_results: bool,
                     ici_mode: str,
                     poly_deg: int,
                     poly_raw: bool,
                     plot_equations: bool,
                     mode: str,
                     all_rows_accum: list,
                     per_run_excel_rows: list,
                     all_eq_rows_accum: list,
                     per_run_eq_rows: list):
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return {"status": "SKIP", "reason": "missing excel"}

    source_label = pretty_source_name(run_root.name)
    tag = tag_from_freq(target_freq)

    if clean_results and mode == "NORM":
        if results_dir.exists():
            backup_results_dir(results_dir)
        results_dir.mkdir(parents=True, exist_ok=True)
    else:
        results_dir.mkdir(parents=True, exist_ok=True)

    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df, sheet_filter=sheet_filter)
    if sip_df.empty:
        return {"status": "FAIL", "reason": "no sip sheets"}

    sat_map, used_sheet, overlap_cnt, block_rows = read_logbook_bestblock_by_overlap(
        log_xlsx, sip_df["measurement_norm"].tolist()
    )
    sip_df["S"] = sip_df["measurement_norm"].map(sat_map).astype(float)

    sip_df["logS"]    = sip_df["S"].apply(safe_log10)
    sip_df["logRho"]  = sip_df["rho"].apply(safe_log10)
    sip_df["logSig"]  = sip_df["sigma_imag"].apply(safe_log10)

    sip_df.copy().to_csv(results_dir / f"debug_before_filter_{tag}Hz__{mode}.csv", index=False, na_rep="")

    mask = np.isfinite(sip_df["logS"]) & np.isfinite(sip_df["logRho"]) & np.isfinite(sip_df["logSig"])
    dropped = sip_df.loc[~mask, ["measurement", "S", "rho", "sigma_imag"]].copy()
    if len(dropped) > 0:
        dropped.to_csv(results_dir / f"dropped_rows_{tag}Hz__{mode}.csv", index=False, na_rep="")
    sip_df = sip_df.loc[mask].copy()

    if len(sip_df) < min_points:
        return {"status": "FAIL", "reason": f"not enough points ({len(sip_df)})"}

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

    rho_slope, rho_int, rho_r2, rho_p, rho_ss, rho_aic = linear_fit(sip_df["logS"], sip_df["logRho"])
    sig_slope, sig_int, sig_r2, sig_p, sig_ss, sig_aic = linear_fit(sip_df["logS"], sip_df["logSig"])

    sip_df["logRho_lin"] = rho_slope*sip_df["logS"] + rho_int
    sip_df["logSig_lin"] = sig_slope*sip_df["logS"] + sig_int

    sip_df["logRho_poly"] = np.nan
    sip_df["logSig_poly"] = np.nan
    r2_logrho_poly = np.nan
    r2_logsig_poly = np.nan
    cR = None
    cS = None
    if poly_deg and poly_deg >= 2:
        cR, _, r2pR, _ = poly_fit_predict(sip_df["logS"], sip_df["logRho"], poly_deg)
        cS, _, r2pS, _ = poly_fit_predict(sip_df["logS"], sip_df["logSig"], poly_deg)
        r2_logrho_poly = r2pR
        r2_logsig_poly = r2pS
        if cR is not None:
            sip_df["logRho_poly"] = np.polyval(cR, sip_df["logS"].to_numpy(float))
        if cS is not None:
            sip_df["logSig_poly"] = np.polyval(cS, sip_df["logS"].to_numpy(float))

    ri_slope, ri_int, ri_r2, ri_p, ri_ss, ri_aic = linear_fit(sip_df["logS"], sip_df["logRhoIdx"])
    sip_df["logRhoIdx_lin"] = ri_slope*sip_df["logS"] + ri_int

    ici_slope, ici_int, ici_r2, ici_p, ici_ss, ici_aic = linear_fit(sip_df["logS"], sip_df["logICI_vn"])
    sip_df["logICI_lin"] = ici_slope*sip_df["logS"] + ici_int

    goal_ri_a, goal_ri_r2, goal_ri_sse = linear_fit_through_origin(sip_df["logS"], sip_df["logRI_vn"])
    goal_ici_a, goal_ici_r2, goal_ici_sse = linear_fit_through_origin(sip_df["logS"], sip_df["logICI_vn"])
    sip_df["logRI_goal"]  = goal_ri_a  * sip_df["logS"]
    sip_df["logICI_goal"] = goal_ici_a * sip_df["logS"]

    eq_raw_rho_lin = eq_linear_str(rho_slope, rho_int, x="logS", y="logRho")
    eq_raw_sig_lin = eq_linear_str(sig_slope, sig_int, x="logS", y="logSig")
    eq_index_ri_lin  = eq_linear_str(ri_slope,  ri_int,  x="logS", y="logRI")
    eq_index_ici_lin = eq_linear_str(ici_slope, ici_int, x="logS", y="logICI")
    eq_target_ri_goal  = eq_origin_str(goal_ri_a,  x="logS", y="logRI")
    eq_target_ici_goal = eq_origin_str(goal_ici_a, x="logS", y="logICI")
    eq_target_raw_n = eq_linear_str(ri_slope,  ri_int,  x="logS", y="logRI")
    eq_target_raw_p = eq_linear_str(ici_slope, ici_int, x="logS", y="logICI")

    eq_logrho_poly = ""
    eq_logsig_poly = ""
    if cR is not None:
        eq_logrho_poly = eq_poly_str(cR, x="logS", y="logRho")
    if cS is not None:
        eq_logsig_poly = eq_poly_str(cS, x="logS", y="logSig")

    sip_df["logRho_log4"]    = np.nan
    sip_df["logRho_logfit"]  = np.nan
    sip_df["resid_log_rho"]  = np.nan

    log4R = None
    has_log_points = False
    log_fit_used = None

    log4_accepted = False
    log4_reason = "not_run"
    eq_logrho_log4 = ""

    if do_log4 and len(sip_df) >= 4:
        log4R, st = logistic_fit_noscipy(
            sip_df["logS"], sip_df["logRho"],
            fast=logistic_fast, max_seconds=logistic_max_seconds
        )
        log4_reason = st

        if log4R is not None:
            ylog4 = logistic4(
                sip_df["logS"].to_numpy(float),
                log4R["L"], log4R["k"], log4R["x0"], log4R["c"]
            )
            sip_df["logRho_log4"]   = ylog4
            sip_df["logRho_logfit"] = ylog4
            sip_df["resid_log_rho"] = sip_df["logRho"] - sip_df["logRho_logfit"]

            eq_logrho_log4 = eq_log4_str(
                log4R["L"], log4R["k"], log4R["x0"], log4R["c"],
                x="logS", y="logRho"
            )

            has_log_points = bool(
                np.isfinite(pd.to_numeric(sip_df["logRho_logfit"], errors="coerce").to_numpy(float)).any()
            )

            ok, why = accept_logistic(rho_r2, rho_aic, log4R)
            log4_accepted = bool(ok)
            log4_reason = why if ok else f"rejected:{why}"
            log_fit_used = log4R

    sip_df["logRho_linfit"] = sip_df["logRho_lin"]
    sip_df["resid_lin_rho"] = sip_df["logRho"] - sip_df["logRho_linfit"]

    sip_df["logSig_log4"] = np.nan
    log4S_log = None
    eq_logsig_log4 = ""
    if do_log4 and len(sip_df) >= 4:
        log4S_log, _ = logistic_fit_noscipy(
            sip_df["logS"], sip_df["logSig"],
            fast=logistic_fast, max_seconds=logistic_max_seconds
        )
        if log4S_log is not None:
            sip_df["logSig_log4"] = logistic4(
                sip_df["logS"].to_numpy(float),
                log4S_log["L"], log4S_log["k"], log4S_log["x0"], log4S_log["c"]
            )
            eq_logsig_log4 = eq_log4_str(
                log4S_log["L"], log4S_log["k"], log4S_log["x0"], log4S_log["c"],
                x="logS", y="logSig"
            )

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

    eq_rho_real_lin = eq_linear_str(mr, cr, x="S", y=("rho_norm" if mode == "NORM" else "rho"))
    eq_sig_real_lin = eq_linear_str(ms, cs, x="S", y=("sig_norm" if mode == "NORM" else "sig"))

    sip_df["rho_poly_real"] = np.nan
    sip_df["sig_poly_real"] = np.nan
    r2_rho_poly_real = np.nan
    r2_sig_poly_real = np.nan
    cpr = None
    cps = None
    eq_rho_real_poly = ""
    eq_sig_real_poly = ""
    if poly_deg and poly_deg >= 2:
        cpr, _, r2pr, _ = poly_fit_predict(xS, y_rho_real, poly_deg)
        cps, _, r2ps, _ = poly_fit_predict(xS, y_sig_real, poly_deg)
        r2_rho_poly_real = r2pr
        r2_sig_poly_real = r2ps
        if cpr is not None:
            sip_df["rho_poly_real"] = np.polyval(cpr, xS)
            eq_rho_real_poly = eq_poly_str(cpr, x="S", y=("rho_norm" if mode == "NORM" else "rho"))
        if cps is not None:
            sip_df["sig_poly_real"] = np.polyval(cps, xS)
            eq_sig_real_poly = eq_poly_str(cps, x="S", y=("sig_norm" if mode == "NORM" else "sig"))

    sip_df["rho_log4_real"] = np.nan
    sip_df["sig_log4_real"] = np.nan
    log4R_real = None
    log4S_real = None
    eq_rho_real_log4 = ""
    eq_sig_real_log4 = ""

    if do_log4 and len(sip_df) >= 4:
        log4R_real, _ = logistic_fit_noscipy(xS, y_rho_real, fast=logistic_fast, max_seconds=logistic_max_seconds)
        log4S_real, _ = logistic_fit_noscipy(xS, y_sig_real, fast=logistic_fast, max_seconds=logistic_max_seconds)
        if log4R_real is not None:
            sip_df["rho_log4_real"] = logistic4(xS, log4R_real["L"], log4R_real["k"], log4R_real["x0"], log4R_real["c"])
            eq_rho_real_log4 = eq_log4_str(
                log4R_real["L"], log4R_real["k"], log4R_real["x0"], log4R_real["c"],
                x="S", y=("rho_norm" if mode == "NORM" else "rho")
            )
        if log4S_real is not None:
            sip_df["sig_log4_real"] = logistic4(xS, log4S_real["L"], log4S_real["k"], log4S_real["x0"], log4S_real["c"])
            eq_sig_real_log4 = eq_log4_str(
                log4S_real["L"], log4S_real["k"], log4S_real["x0"], log4S_real["c"],
                x="S", y=("sig_norm" if mode == "NORM" else "sig")
            )

    joined_name = f"joined__{mode}_{tag}Hz.csv"
    joined_path = results_dir / joined_name
    sip_df.to_csv(joined_path, index=False, na_rep="")

    if mode == "NORM":
        (results_dir / f"joined_{tag}Hz.csv").write_text(joined_path.read_text())

    if mode == "NORM":
        legacy_csv = f"joined_{tag}Hz.csv"
        write_plot_raw_exponents(
            results_dir, legacy_csv, f"exponents_raw_{tag}Hz.png", target_freq,
            source_label,
            rho_slope, rho_int, rho_r2, rho_p, eq_raw_rho_lin,
            sig_slope, sig_int, sig_r2, sig_p, eq_raw_sig_lin,
            show_eq=plot_equations
        )
        write_plot_index_exponents(
            results_dir, legacy_csv, f"exponents_index_{tag}Hz.png", target_freq,
            source_label,
            ri_slope, ri_int, ri_r2, ri_p, eq_index_ri_lin,
            ici_slope, ici_int, ici_r2, ici_p, eq_index_ici_lin,
            show_eq=plot_equations
        )
        write_plot_target_exponents(
            results_dir, legacy_csv, f"target_exponents_{tag}Hz.png", target_freq,
            source_label,
            goal_ri_a, goal_ri_r2, goal_ri_sse, eq_target_ri_goal,
            goal_ici_a, goal_ici_r2, goal_ici_sse, eq_target_ici_goal,
            show_eq=plot_equations
        )
        write_plot_target_raw_exponents(
            results_dir, legacy_csv, f"target_raw_exponents_{tag}Hz.png", target_freq,
            source_label,
            -ri_slope, ri_r2, eq_target_raw_n,
            -ici_slope, ici_r2, eq_target_raw_p,
            show_eq=plot_equations
        )
        write_plot_logistic_diag(
            results_dir, legacy_csv, f"exponents_logistic_diagnostics_{tag}Hz.png",
            target_freq, source_label,
            rho_slope, rho_int, rho_r2, eq_raw_rho_lin,
            has_log_points=has_log_points,
            log4_params=log_fit_used,
            log4_eq=eq_logrho_log4,
            log4_accepted=log4_accepted,
            log4_reason=log4_reason,
            show_eq=plot_equations
        )
        write_plot_logistic_only(
            results_dir, legacy_csv, f"exponents_logistic_only_{tag}Hz.png",
            target_freq, source_label,
            has_log_points=has_log_points,
            log4_eq=eq_logrho_log4,
            log4_accepted=log4_accepted,
            log4_reason=log4_reason,
            show_eq=plot_equations
        )

        for sc in [
            "plot_raw_exponents.gp",
            "plot_index_exponents.gp",
            "plot_target_exponents.gp",
            "plot_target_raw_exponents.gp",
            "plot_logistic_diagnostics.gp",
            "plot_logistic_only.gp"
        ]:
            _run_gnuplot(sc, cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

    methods = "LINEAR"
    if poly_deg and poly_deg >= 2:
        methods += "+POLY"
    if do_log4:
        methods += "+LOG4"

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
        y_log4_col="logRho_log4",
        linear_eq=eq_raw_rho_lin,
        poly_eq=eq_logrho_poly,
        log4_eq=eq_logrho_log4,
        show_eq=plot_equations
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
        y_log4_col="logSig_log4",
        linear_eq=eq_raw_sig_lin,
        poly_eq=eq_logsig_poly,
        log4_eq=eq_logsig_log4,
        show_eq=plot_equations
    )
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
        y_log4_col="rho_log4_real",
        linear_eq=eq_rho_real_lin,
        poly_eq=eq_rho_real_poly,
        log4_eq=eq_rho_real_log4,
        show_eq=plot_equations
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
        y_log4_col="sig_log4_real",
        linear_eq=eq_sig_real_lin,
        poly_eq=eq_sig_real_poly,
        log4_eq=eq_sig_real_log4,
        show_eq=plot_equations
    )

    for sc in [
        f"rho_vs_S__LOGLOG__{mode}__{methods}_{tag}Hz.gp",
        f"cond_vs_S__LOGLOG__{mode}__{methods}_{tag}Hz.gp",
        f"rho_vs_S__REAL__{mode}__{methods}_{tag}Hz.gp",
        f"cond_vs_S__REAL__{mode}__{methods}_{tag}Hz.gp",
    ]:
        _run_gnuplot(sc, cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)

    for eq_list in (all_eq_rows_accum, per_run_eq_rows):
        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "exponents_raw", "raw_rho_linear", eq_raw_rho_lin, rho_r2, rho_p, rho_aic)
        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "exponents_raw", "raw_sig_linear", eq_raw_sig_lin, sig_r2, sig_p, sig_aic)

        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "exponents_index", "index_ri_linear", eq_index_ri_lin, ri_r2, ri_p, ri_aic)
        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "exponents_index", "index_ici_linear", eq_index_ici_lin, ici_r2, ici_p, ici_aic)

        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "target_exponents", "target_ri_origin", eq_target_ri_goal, goal_ri_r2, np.nan, np.nan)
        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "target_exponents", "target_ici_origin", eq_target_ici_goal, goal_ici_r2, np.nan, np.nan)

        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "target_raw_exponents", "target_raw_n_linear", eq_target_raw_n, ri_r2, ri_p, ri_aic)
        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "target_raw_exponents", "target_raw_p_linear", eq_target_raw_p, ici_r2, ici_p, ici_aic)

        if eq_logrho_poly:
            add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                             "rho_vs_S__LOGLOG", f"logrho_poly_deg{poly_deg}", eq_logrho_poly, r2_logrho_poly)
        if eq_logsig_poly:
            add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                             "cond_vs_S__LOGLOG", f"logsig_poly_deg{poly_deg}", eq_logsig_poly, r2_logsig_poly)
        if eq_logrho_log4:
            add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                             "rho_vs_S__LOGLOG", "logrho_logistic4", eq_logrho_log4,
                             log4R.get("r2", np.nan) if log4R else np.nan,
                             log4R.get("p_model", np.nan) if log4R else np.nan,
                             log4R.get("aic", np.nan) if log4R else np.nan,
                             "accepted" if log4_accepted else "rejected",
                             log4_reason)
        if eq_logsig_log4:
            add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                             "cond_vs_S__LOGLOG", "logsig_logistic4", eq_logsig_log4,
                             log4S_log.get("r2", np.nan) if log4S_log else np.nan,
                             log4S_log.get("p_model", np.nan) if log4S_log else np.nan,
                             log4S_log.get("aic", np.nan) if log4S_log else np.nan)

        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "rho_vs_S__REAL", "rho_real_linear", eq_rho_real_lin, r2r, pr, aicr)
        add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                         "cond_vs_S__REAL", "sig_real_linear", eq_sig_real_lin, r2s, ps, aics)

        if eq_rho_real_poly:
            add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                             "rho_vs_S__REAL", f"rho_real_poly_deg{poly_deg}", eq_rho_real_poly, r2_rho_poly_real)
        if eq_sig_real_poly:
            add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                             "cond_vs_S__REAL", f"sig_real_poly_deg{poly_deg}", eq_sig_real_poly, r2_sig_poly_real)
        if eq_rho_real_log4:
            add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                             "rho_vs_S__REAL", "rho_real_logistic4", eq_rho_real_log4,
                             log4R_real.get("r2", np.nan) if log4R_real else np.nan,
                             log4R_real.get("p_model", np.nan) if log4R_real else np.nan,
                             log4R_real.get("aic", np.nan) if log4R_real else np.nan)
        if eq_sig_real_log4:
            add_equation_row(eq_list, run_root, source_label, mode, target_freq,
                             "cond_vs_S__REAL", "sig_real_logistic4", eq_sig_real_log4,
                             log4S_real.get("r2", np.nan) if log4S_real else np.nan,
                             log4S_real.get("p_model", np.nan) if log4S_real else np.nan,
                             log4S_real.get("aic", np.nan) if log4S_real else np.nan)

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
        row["log4_diag_status"] = ("accepted" if log4_accepted else "rejected") if (log_fit_used is not None) else "no_fit"
        row["log4_diag_reason"] = str(log4_reason)
        row["eq_raw_rho_lin"] = eq_raw_rho_lin
        row["eq_raw_sig_lin"] = eq_raw_sig_lin
        row["eq_index_ri_lin"] = eq_index_ri_lin
        row["eq_index_ici_lin"] = eq_index_ici_lin
        row["eq_target_ri_goal"] = eq_target_ri_goal
        row["eq_target_ici_goal"] = eq_target_ici_goal
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
        "plot_equations": bool(plot_equations),
        "has_log_points": bool(has_log_points),
        "log4_diag_status": ("accepted" if log4_accepted else "rejected") if (log_fit_used is not None) else "no_fit",
        "log4_diag_reason": str(log4_reason),
        "r2_raw_rho_lin": float(rho_r2) if np.isfinite(rho_r2) else rho_r2,
        "r2_raw_sig_lin": float(sig_r2) if np.isfinite(sig_r2) else sig_r2,
        "r2_index_ri_lin": float(ri_r2) if np.isfinite(ri_r2) else ri_r2,
        "r2_index_ici_lin": float(ici_r2) if np.isfinite(ici_r2) else ici_r2,
        "goal_ri_a": float(goal_ri_a) if np.isfinite(goal_ri_a) else goal_ri_a,
        "goal_ici_a": float(goal_ici_a) if np.isfinite(goal_ici_a) else goal_ici_a,
        "eq_raw_rho_lin": eq_raw_rho_lin,
        "eq_raw_sig_lin": eq_raw_sig_lin,
        "eq_index_ri_lin": eq_index_ri_lin,
        "eq_index_ici_lin": eq_index_ici_lin,
        "eq_target_ri_goal": eq_target_ri_goal,
        "eq_target_ici_goal": eq_target_ici_goal,
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

    ap.add_argument("--clean-results", action="store_true",
                    help="BACKUP existing Results/ to Results__backup_<timestamp>/ before writing new outputs.")

    ap.add_argument("--ici-mode", type=str, default="sigref_over_sig",
                    choices=["sigref_over_sig","sig_over_sigref"],
                    help="Affects normalized index definition (NORM mode). VN TARGET always uses sigref/sig.")

    ap.add_argument("--poly-deg", type=int, default=0, help="Polynomial degree (>=2). 0 disables.")
    ap.add_argument("--poly-raw", action="store_true", help="Kept for compatibility (poly already applied).")

    ap.add_argument("--logistic-fast", action="store_true", help="Use SciPy-free fast logistic search")
    ap.add_argument("--logistic-max-seconds", type=float, default=8.0, help="Time cap for logistic per run")
    ap.add_argument("--logistic4", action="store_true", help="Enable 4-parameter logistic fits (LOG4).")

    ap.add_argument("--plot-equations", action="store_true",
                    help="Write fit equations on plots. Default is OFF. Excel equations sheet is always written.")

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
    all_eq_rows = []
    summary_rows = []
    t0 = time.time()

    for i, run in enumerate(roots, 1):
        print(f"[{i}/{len(roots)}] START {run}", flush=True)
        ti = time.time()

        per_run_excel_rows = []
        per_run_eq_rows = []

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
                    plot_equations=args.plot_equations,
                    mode=mode,
                    all_rows_accum=all_rows,
                    per_run_excel_rows=per_run_excel_rows,
                    all_eq_rows_accum=all_eq_rows,
                    per_run_eq_rows=per_run_eq_rows
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
                df_eq  = pd.DataFrame(per_run_eq_rows)
                out_xlsx = (run / RESULTS_DIR) / f"all_values_{tag}Hz.xlsx"
                with pd.ExcelWriter(out_xlsx, engine="openpyxl") as w:
                    df_run.to_excel(w, sheet_name="all_values", index=False)
                    if not df_eq.empty:
                        df_eq.to_excel(w, sheet_name="equations", index=False)

        except subprocess.TimeoutExpired:
            summary_rows.append({"run": str(run), "status": "FAIL", "reason": "gnuplot timeout", "points": 0})
            print(f"[{i}/{len(roots)}] FAIL gnuplot timeout", flush=True)
        except Exception as e:
            summary_rows.append({"run": str(run), "status": "FAIL", "reason": str(e), "points": 0})
            print(f"[{i}/{len(roots)}] FAIL {e}", flush=True)

    if summary_rows:
        df_sum = pd.DataFrame(summary_rows)
        sum_csv = base / f"Results_Summary_{tag}Hz.csv"
        sum_xlsx = base / f"Results_Summary_{tag}Hz.xlsx"
        df_sum.to_csv(sum_csv, index=False)
        with pd.ExcelWriter(sum_xlsx, engine="openpyxl") as w:
            df_sum.to_excel(w, sheet_name="summary", index=False)
            if all_eq_rows:
                pd.DataFrame(all_eq_rows).to_excel(w, sheet_name="equations", index=False)
        print(f"\nSummary written: {sum_csv}", flush=True)
        print(f"Summary written: {sum_xlsx}", flush=True)

    if all_rows:
        df_all = pd.DataFrame(all_rows)
        df_eq_all = pd.DataFrame(all_eq_rows)
        all_xlsx = base / f"all_values_{tag}Hz.xlsx"
        with pd.ExcelWriter(all_xlsx, engine="openpyxl") as w:
            df_all.to_excel(w, sheet_name="all_values", index=False)
            if not df_eq_all.empty:
                df_eq_all.to_excel(w, sheet_name="equations", index=False)
        print(f"All-values written: {all_xlsx}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)

if __name__ == "__main__":
    main()

# =============================================================================
# RUNNER
# =============================================================================
# Default (equations only in Excel, NOT on plots):
# python3 sip_pipeline.py "/mnt/3rd900/Projects/LA Project new/For_LBNL/SIP" \
#   --freq 0.01 --min-points 4 \
#   --logistic-fast --logistic-max-seconds 8 \
#   --clean-results --both-modes \
#   --ici-mode sig_over_sigref \
#   --poly-deg 2 --poly-raw \
#   --logistic4
#
# With equations on plots too:
# python3 sip_pipeline.py "/mnt/3rd900/Projects/LA Project new/For_LBNL/SIP" \
#   --freq 0.01 --min-points 4 \
#   --logistic-fast --logistic-max-seconds 8 \
#   --clean-results --both-modes \
#   --ici-mode sig_over_sigref \
#   --poly-deg 2 --poly-raw \
#   --logistic4 \
#   --plot-equations
