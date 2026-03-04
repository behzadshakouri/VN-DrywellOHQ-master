#!/usr/bin/env python3
"""
auto_sip_pipeline_v27_noscipy.py  (v27+++)

KEY PROMISE (THIS VERSION)
--------------------------
✅ Keeps the ORIGINAL v27 outputs (unchanged names/logic) for the *normalized* run:
  Results/joined_<tag>Hz.csv
  Results/fit_report_<tag>Hz.txt
  Results/target_exponents_<tag>Hz.png
  Results/index_exponents_<tag>Hz.png
  Results/exponents_raw_<tag>Hz.png
  Results/exponents_logistic_only_<tag>Hz.png
  Results/exponents_logistic_diagnostics_<tag>Hz.png
  Results/debug_before_filter_<tag>Hz.csv
  Results/dropped_rows_<tag>Hz.csv (if any)

✅ Adds NEW capabilities, WITHOUT breaking legacy:
  - run BOTH modes: NORM and UNORM (use --both-modes)
  - run BOTH spaces for NEW analyses: REAL-REAL and LOG-LOG
  - run ALL methods for NEW analyses: linear, polynomial (deg>=2), logistic4 (4-param) in both spaces and both modes
  - write per-run Excel with all values/logs/calcs/fits + global Excel summary
  - add "no-normalization" (UNORM) versions of plots/files using clear suffixes, so no overwrites

FINAL GOAL SUPPORT
------------------
The central new set of fits/plots is conductivity vs moisture:
  sigma_imag vs S   (REAL and LOGLOG; NORM and UNORM; linear/poly/logistic4)

Normalization meaning:
  - NORM: also compute RI/ICI indices (as before) and VN target (as before)
  - UNORM: NO reference normalization; focus on direct rho, sigma_imag relations vs S

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

# Logistic bounds / acceptance (legacy)
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
    s = re.sub(r'\bS-(\d+)\b', r'S_{-\1}', s)
    return s

def tag_freq(freq: float) -> str:
    return str(freq).replace(".", "p")

def mode_suffix(normalize: bool) -> str:
    return "NORM" if normalize else "UNORM"

def space_suffix(is_log: bool) -> str:
    return "LOGLOG" if is_log else "REAL"

def fname(prefix: str, normalize: bool, is_log: bool, tag: str, extra: str = "", ext: str = ".png") -> str:
    base = f"{prefix}__{space_suffix(is_log)}__{mode_suffix(normalize)}"
    if extra:
        base += f"__{extra}"
    base += f"_{tag}Hz{ext}"
    return base


# =============================================================================
# F-DISTRIBUTION + P-VALUE (NO SCIPY)  (legacy)
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
# FIT ROUTINES (legacy + used in new)
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
# 4-PARAM LOGISTIC FIT (legacy logistic4 + now applied broadly)
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
        "ss_res": float(ss_res), "ss_tot": float(ss_tot),
        "n": int(n)
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
    if not results_dir.exists():
        return
    for p in results_dir.iterdir():
        try:
            if p.is_file() or p.is_symlink():
                p.unlink()
        except Exception:
            pass


# =============================================================================
# LEGACY GNUPLOT WRITERS (UNCHANGED NAMES/STYLE)
# =============================================================================

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

def write_plot_logistic_only(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                             source_label: str, log_fit_acc):
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
                             source_label: str, has_log: bool):
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
# NEW: SIMPLE XY PLOT WRITER (for new REAL/LOGLOG and UNORM outputs)
# =============================================================================

def write_plot_xy(results_dir: Path,
                  csv_name: str,
                  out_png: str,
                  source_label: str,
                  title: str,
                  xlabel: str,
                  ylabel: str,
                  x_col: str,
                  y_col: str,
                  y_lin_col: str,
                  y_poly_col: str|None,
                  y_log4_col: str|None):
    poly_clause = ""
    if y_poly_col:
        poly_clause = f", \\\n     '{csv_name}' using (column(\"{x_col}\")):(column(\"{y_poly_col}\")) with lines dt 1 lw 4 title 'poly'"
    log_clause = ""
    if y_log4_col:
        log_clause = f", \\\n     '{csv_name}' using (column(\"{x_col}\")):(column(\"{y_log4_col}\")) with lines dt 1 lw 4 title 'logistic4'"

    gp_name = Path(out_png).stem + ".gp"
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

plot '{csv_name}' using (column("{x_col}")):(column("{y_col}")) with points pt 7 ps 1.8 title 'data', \\
     '{csv_name}' using (column("{x_col}")):(column("{y_lin_col}")) with lines dt 2 lw 4 title 'linear'{poly_clause}{log_clause}
"""
    (results_dir / gp_name).write_text(gp)


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
# MODE APPLICATION (NORM vs UNORM)
# =============================================================================

def apply_mode(df_in: pd.DataFrame, normalize: bool, ici_mode: str):
    """
    Build derived columns for the selected mode.
    Always creates:
      logS, logRho, logImag
    If normalize:
      compute references and GOAL/INDEX + VN target columns (legacy definitions)
    """
    df = df_in.copy()
    df["logS"]    = df["S"].apply(safe_log10)
    df["logRho"]  = df["rho"].apply(safe_log10)
    df["logImag"] = df["sigma_imag"].apply(safe_log10)

    if normalize:
        iref = int(np.nanargmax(df["S"].to_numpy(float)))
        rho_ref  = float(df.iloc[iref]["rho"])
        sig_ref  = float(df.iloc[iref]["sigma_imag"])
        logRho_ref = safe_log10(rho_ref)
        logSig_ref = safe_log10(sig_ref)

        df["rho_ref"] = rho_ref
        df["sig_ref"] = sig_ref
        df["logRho_ref"] = logRho_ref
        df["logSig_ref"] = logSig_ref

        df["rho_idx"]   = df["rho"] / rho_ref
        df["logRhoIdx"] = df["logRho"] - logRho_ref

        if ici_mode == "sig_over_sigref":
            df["imag_idx"]   = df["sigma_imag"] / sig_ref
            df["logImagIdx"] = df["logImag"] - logSig_ref
        else:
            df["imag_idx"]   = sig_ref / df["sigma_imag"]
            df["logImagIdx"] = logSig_ref - df["logImag"]

        # VN target (always downtrend)
        df["logRI_vn"]  = df["logRhoIdx"]
        df["logICI_vn"] = logSig_ref - df["logImag"]
    else:
        # UNORM: still has logs of rho and sigma (no indices)
        pass

    return df


# =============================================================================
# FITTING HELPERS FOR NEW FEATURES
# =============================================================================

def fit_all_methods(x, y, poly_deg: int, do_log4: bool,
                    logistic_min_points: int, logistic_fast: bool, logistic_max_seconds: float):
    """
    Fit linear + optional poly + optional logistic4 (accepted vs linear).
    Returns dict with:
      lin: slope, intercept, r2, p, aic
      poly: coeffs, r2 (or None)
      log4_acc: dict or None
      log4_status, log4_reason
    """
    slope, intercept, r2, p, _, aic = linear_fit(x, y)
    lin = {"slope": slope, "intercept": intercept, "r2": r2, "p": p, "aic": aic}

    poly = {"coeffs": None, "r2": np.nan}
    if poly_deg >= 2:
        coeffs, _, r2p, _ = poly_fit_predict(x, y, poly_deg)
        poly = {"coeffs": coeffs, "r2": r2p}

    log4 = None
    log4_status = "skipped"
    log4_reason = "n/a"
    log4_acc = None
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

def add_fit_cols(df: pd.DataFrame, x_col: str, y_col: str,
                 y_lin: str, y_poly: str|None, y_log4: str|None,
                 lin: dict, poly: dict, log4_acc: dict|None):
    m = np.isfinite(df[x_col]) & np.isfinite(df[y_col])
    xx = df.loc[m, x_col].to_numpy(float)
    df.loc[m, y_lin] = lin["slope"]*xx + lin["intercept"]
    if y_poly and (poly.get("coeffs") is not None):
        df.loc[m, y_poly] = np.polyval(poly["coeffs"], xx)
    if y_log4 and (log4_acc is not None):
        L = log4_acc["L"]; k = log4_acc["k"]; x0 = log4_acc["x0"]; c = log4_acc["c"]
        df.loc[m, y_log4] = logistic4(xx, L, k, x0, c)
    return df


# =============================================================================
# MAIN PROCESSING FOR ONE RUN (LEGACY + NEW)
# =============================================================================

def process_run(run_root: Path, target_freq: float, sheet_filter: str|None,
                dry_run: bool, min_points: int, max_df: float,
                do_logistic_legacy: bool, logistic_min_points: int,
                logistic_fast: bool, logistic_max_seconds: float,
                gnuplot_timeout: int, clean_results: bool,
                ici_mode: str, both_modes: bool,
                poly_deg: int, poly_raw: bool,
                do_logistic4: bool):
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

    tag = tag_freq(target_freq)

    # Read SIP sheets
    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df, sheet_filter=sheet_filter)
    if sip_df.empty:
        return [{
            "run": str(run_root), "source": source_label, "mode": "N/A",
            "status": "FAIL", "reason": "No SIP sheets readable / target freq missing", "points": 0
        }]

    # Read logbook and map saturation
    sat_map, used_sheet, overlap_cnt, block_rows = read_logbook_bestblock_by_overlap(
        log_xlsx, sip_df["measurement_norm"].tolist()
    )
    sip_df["S"] = sip_df["measurement_norm"].map(sat_map).astype(float)

    # Compute base logs (for debug) and write debug-before-filter (legacy name)
    tmp_dbg = sip_df.copy()
    tmp_dbg["logS"] = tmp_dbg["S"].apply(safe_log10)
    tmp_dbg["logRho"] = tmp_dbg["rho"].apply(safe_log10)
    tmp_dbg["logImag"] = tmp_dbg["sigma_imag"].apply(safe_log10)
    tmp_dbg.to_csv(results_dir / f"debug_before_filter_{tag}Hz.csv", index=False)

    # Filter invalid (legacy filter: require finite logS/logRho/logImag)
    mask = np.isfinite(tmp_dbg["logS"]) & np.isfinite(tmp_dbg["logRho"]) & np.isfinite(tmp_dbg["logImag"])
    dropped = tmp_dbg.loc[~mask, ["measurement", "S", "rho", "sigma_imag"]].copy()
    if len(dropped) > 0:
        dropped.to_csv(results_dir / f"dropped_rows_{tag}Hz.csv", index=False)
    sip_df = sip_df.loc[mask].copy()

    if len(sip_df) < min_points:
        return [{
            "run": str(run_root), "source": source_label, "mode": "N/A",
            "status": "FAIL",
            "reason": f"Not enough matched points (need >= {min_points})",
            "points": int(len(sip_df)),
            "logbook_sheet": used_sheet,
            "logbook_overlap": overlap_cnt,
            "logbook_block_rows": block_rows,
        }]

    # We will run modes:
    #   normalize=True  => legacy v27 behavior (and keep original filenames)
    #   normalize=False => UNORM extra outputs with suffixes (no overwrite)
    modes = [True] if not both_modes else [True, False]

    # Per-run excel aggregation
    excel_path = results_dir / f"all_values_{tag}Hz.xlsx"
    excel_sheets = {}
    fit_rows = []

    out_rows = []

    for normalize in modes:
        df = apply_mode(sip_df, normalize=normalize, ici_mode=ici_mode)

        # -----------------------------
        # LEGACY v27 pipeline (ONLY for normalize=True)
        # -----------------------------
        if normalize:
            # reference point as max S (legacy)
            iref = int(np.nanargmax(df["S"].to_numpy(float)))
            rho_ref  = float(df.iloc[iref]["rho"])
            sig_ref  = float(df.iloc[iref]["sigma_imag"])
            logRho_ref = safe_log10(rho_ref)
            logSig_ref = safe_log10(sig_ref)

            # RAW linear fits (legacy)
            slope_r, intercept_r, r2_r, p_r, ss_res_lin_r, aic_lin_r = linear_fit(df["logS"], df["logRho"])
            slope_i, intercept_i, r2_i, p_i, ss_res_lin_i, aic_lin_i = linear_fit(df["logS"], df["logImag"])
            lin_rho = {"slope": slope_r, "intercept": intercept_r, "r2": r2_r, "p": p_r, "aic": aic_lin_r}
            lin_im  = {"slope": slope_i, "intercept": intercept_i, "r2": r2_i, "p": p_i, "aic": aic_lin_i}

            # GOAL/INDEX through-origin (legacy)
            a_ri,  r2_ri,  _ = linear_fit_through_origin(df["logS"], df["logRhoIdx"])
            a_ici, r2_ici, _ = linear_fit_through_origin(df["logS"], df["logImagIdx"])

            n_idx = -a_ri
            if ici_mode == "sig_over_sigref":
                p_idx = +a_ici
            else:
                p_idx = -a_ici

            # TARGET/VN with intercept (legacy)
            slope_RI_vn,  int_RI_vn,  r2_RI_vn,  _, _, _ = linear_fit(df["logS"], df["logRI_vn"])
            slope_ICI_vn, int_ICI_vn, r2_ICI_vn, _, _, _ = linear_fit(df["logS"], df["logICI_vn"])
            n_vn = -slope_RI_vn
            p_vn = -slope_ICI_vn

            df["logRI_vn_fit"]  = slope_RI_vn  * df["logS"] + int_RI_vn
            df["logICI_vn_fit"] = slope_ICI_vn * df["logS"] + int_ICI_vn

            # LEGACY logistic on RAW logRho vs logS
            log_fit = None
            log_status = "skipped"
            if do_logistic_legacy and len(df) >= logistic_min_points:
                log_fit, log_status = logistic_fit_noscipy(
                    df["logS"], df["logRho"],
                    fast=logistic_fast,
                    max_seconds=(logistic_max_seconds if logistic_max_seconds and logistic_max_seconds > 0 else 0.0)
                )
            elif do_logistic_legacy:
                log_status = "not_enough_points"

            log_fit_acc = None
            log_accept_reason = "n/a"
            if log_fit is not None:
                ok, reason = accept_logistic(r2_r, aic_lin_r, log_fit)
                log_accept_reason = reason
                if ok:
                    log_fit_acc = log_fit

            # Residuals and fitted curves for diagnostics (legacy)
            logS_arr = df["logS"].to_numpy(float)
            logR_arr = df["logRho"].to_numpy(float)
            lin_pred = slope_r*logS_arr + intercept_r
            df["resid_lin_rho"] = logR_arr - lin_pred

            if log_fit_acc is not None:
                L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
                log_pred = logistic4(logS_arr, L, k, x0, c)
                df["resid_log_rho"] = logR_arr - log_pred
                df["logRho_linfit"] = lin_pred
                df["logRho_logfit"] = log_pred
            else:
                df["resid_log_rho"] = np.nan
                df["logRho_linfit"] = lin_pred
                df["logRho_logfit"] = np.nan

            # Write LEGACY joined CSV (exact name)
            csv_path_legacy = results_dir / f"joined_{tag}Hz.csv"
            df_cols_legacy = [
                "measurement","measurement_norm","freq_hz","S","rho","sigma_imag",
                "logS","logRho","logImag",
                "rho_idx","imag_idx","logRhoIdx","logImagIdx",
                "logRI_vn","logICI_vn","logRI_vn_fit","logICI_vn_fit",
                "resid_lin_rho","resid_log_rho",
                "logRho_linfit","logRho_logfit",
            ]
            for c in df_cols_legacy:
                if c not in df.columns:
                    df[c] = np.nan
            df[df_cols_legacy].to_csv(csv_path_legacy, index=False)

            # Write LEGACY report (exact name)
            report = results_dir / f"fit_report_{tag}Hz.txt"
            report.write_text(
                f"Run root: {run_root}\n"
                f"Source (folder): {source_label}\n"
                f"SIP file: {sip_xlsx}\n"
                f"Logbook file: {log_xlsx}\n"
                f"Logbook sheet used: {used_sheet}\n"
                f"Logbook block rows (usable): {block_rows}\n"
                f"Overlap with SIP measurement IDs: {overlap_cnt}\n\n"
                f"Matched points used: {len(df)}\n"
                f"Reference point (max S): {df.iloc[iref]['measurement']}  S_ref={df.iloc[iref]['S']:.6g}\n"
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

            # Write LEGACY plots (exact names)
            write_plot_target(results_dir, csv_path_legacy.name, f"target_exponents_{tag}Hz.png", target_freq,
                              source_label, n_vn, r2_RI_vn, p_vn, r2_ICI_vn)
            write_plot_index(results_dir, csv_path_legacy.name, f"index_exponents_{tag}Hz.png", target_freq,
                             source_label, n_idx, r2_ri, p_idx, r2_ici)
            write_plot_raw(results_dir, csv_path_legacy.name, f"exponents_raw_{tag}Hz.png", target_freq,
                           source_label, lin_rho, lin_im, log_fit_acc)
            write_plot_logistic_only(results_dir, csv_path_legacy.name, f"exponents_logistic_only_{tag}Hz.png",
                                     target_freq, source_label, log_fit_acc)
            write_plot_logistic_diag(results_dir, csv_path_legacy.name, f"exponents_logistic_diagnostics_{tag}Hz.png",
                                     target_freq, source_label, has_log=(log_fit_acc is not None))

            # Execute LEGACY gnuplot unless dry-run
            if not dry_run:
                _run_gnuplot("plot_target_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
                _run_gnuplot("plot_index_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
                _run_gnuplot("plot_raw_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
                _run_gnuplot("plot_logistic_only.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
                _run_gnuplot("plot_logistic_diagnostics.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

        # -----------------------------
        # NEW additions (for BOTH modes): REAL+LOGLOG, linear+poly+logistic4
        # Focus on sigma_imag vs S (goal), plus rho vs S
        # -----------------------------

        # Ensure positivity for REAL fits
        df["S_real"] = df["S"].astype(float)
        df["rho_real"] = df["rho"].astype(float)
        df["sig_real"] = df["sigma_imag"].astype(float)

        # Prepare new fit columns
        for col in [
            "rho_lin_real", "rho_poly_real", "rho_log4_real",
            "sig_lin_real", "sig_poly_real", "sig_log4_real",
            "logRho_lin", "logRho_poly", "logRho_log4",
            "logSig_lin", "logSig_poly", "logSig_log4"
        ]:
            if col not in df.columns:
                df[col] = np.nan

        # REAL fits: rho vs S, sig vs S
        mR = np.isfinite(df["S_real"]) & np.isfinite(df["rho_real"]) & np.isfinite(df["sig_real"]) & (df["S_real"] > 0) & (df["rho_real"] > 0) & (df["sig_real"] > 0)
        dR = df.loc[mR].copy()

        if len(dR) >= min_points:
            # rho vs S
            lin, poly, log4_acc, log4_status, log4_reason = fit_all_methods(
                dR["S_real"], dR["rho_real"], poly_deg=poly_deg,
                do_log4=do_logistic4,
                logistic_min_points=logistic_min_points,
                logistic_fast=logistic_fast,
                logistic_max_seconds=logistic_max_seconds
            )
            df = add_fit_cols(df, "S_real", "rho_real",
                              "rho_lin_real",
                              "rho_poly_real" if poly_deg >= 2 else None,
                              "rho_log4_real" if do_logistic4 else None,
                              lin, poly, log4_acc)

            fit_rows.append({
                "run": str(run_root), "source": source_label, "mode": mode_suffix(normalize),
                "space": "REAL", "x": "S", "y": "rho",
                "lin_slope": lin["slope"], "lin_intercept": lin["intercept"], "lin_r2": lin["r2"], "lin_p": lin["p"], "lin_aic": lin["aic"],
                "poly_deg": poly_deg if poly_deg >= 2 else 0, "poly_r2": poly["r2"],
                "log4_attempted": bool(do_logistic4), "log4_accepted": bool(log4_acc is not None),
                "log4_status": log4_status, "log4_reason": log4_reason,
                "log4_L": (log4_acc["L"] if log4_acc else np.nan),
                "log4_k": (log4_acc["k"] if log4_acc else np.nan),
                "log4_x0": (log4_acc["x0"] if log4_acc else np.nan),
                "log4_c": (log4_acc["c"] if log4_acc else np.nan),
                "log4_r2": (log4_acc["r2"] if log4_acc else np.nan),
                "log4_p": (log4_acc["p_model"] if log4_acc else np.nan),
                "log4_aic": (log4_acc["aic"] if log4_acc else np.nan),
            })

            # sigma vs S (MAIN GOAL)
            linS, polyS, log4S_acc, log4S_status, log4S_reason = fit_all_methods(
                dR["S_real"], dR["sig_real"], poly_deg=poly_deg,
                do_log4=do_logistic4,
                logistic_min_points=logistic_min_points,
                logistic_fast=logistic_fast,
                logistic_max_seconds=logistic_max_seconds
            )
            df = add_fit_cols(df, "S_real", "sig_real",
                              "sig_lin_real",
                              "sig_poly_real" if poly_deg >= 2 else None,
                              "sig_log4_real" if do_logistic4 else None,
                              linS, polyS, log4S_acc)

            fit_rows.append({
                "run": str(run_root), "source": source_label, "mode": mode_suffix(normalize),
                "space": "REAL", "x": "S", "y": "sigma_imag",
                "lin_slope": linS["slope"], "lin_intercept": linS["intercept"], "lin_r2": linS["r2"], "lin_p": linS["p"], "lin_aic": linS["aic"],
                "poly_deg": poly_deg if poly_deg >= 2 else 0, "poly_r2": polyS["r2"],
                "log4_attempted": bool(do_logistic4), "log4_accepted": bool(log4S_acc is not None),
                "log4_status": log4S_status, "log4_reason": log4S_reason,
                "log4_L": (log4S_acc["L"] if log4S_acc else np.nan),
                "log4_k": (log4S_acc["k"] if log4S_acc else np.nan),
                "log4_x0": (log4S_acc["x0"] if log4S_acc else np.nan),
                "log4_c": (log4S_acc["c"] if log4S_acc else np.nan),
                "log4_r2": (log4S_acc["r2"] if log4S_acc else np.nan),
                "log4_p": (log4S_acc["p_model"] if log4S_acc else np.nan),
                "log4_aic": (log4S_acc["aic"] if log4S_acc else np.nan),
            })

        # LOGLOG fits: logRho vs logS, logSig vs logS
        df["logSig"] = df["logImag"]  # naming convenience for new fits
        mL = np.isfinite(df["logS"]) & np.isfinite(df["logRho"]) & np.isfinite(df["logSig"])
        dL = df.loc[mL].copy()

        if len(dL) >= min_points:
            # logRho vs logS
            linLR, polyLR, log4LR_acc, log4LR_status, log4LR_reason = fit_all_methods(
                dL["logS"], dL["logRho"],
                poly_deg=(poly_deg if poly_raw else 0),
                do_log4=do_logistic4,
                logistic_min_points=logistic_min_points,
                logistic_fast=logistic_fast,
                logistic_max_seconds=logistic_max_seconds
            )
            df = add_fit_cols(df, "logS", "logRho",
                              "logRho_lin",
                              "logRho_poly" if (poly_raw and poly_deg >= 2) else None,
                              "logRho_log4" if do_logistic4 else None,
                              linLR, polyLR, log4LR_acc)

            # logSig vs logS
            linLS, polyLS, log4LS_acc, log4LS_status, log4LS_reason = fit_all_methods(
                dL["logS"], dL["logSig"],
                poly_deg=poly_deg,
                do_log4=do_logistic4,
                logistic_min_points=logistic_min_points,
                logistic_fast=logistic_fast,
                logistic_max_seconds=logistic_max_seconds
            )
            df = add_fit_cols(df, "logS", "logSig",
                              "logSig_lin",
                              "logSig_poly" if poly_deg >= 2 else None,
                              "logSig_log4" if do_logistic4 else None,
                              linLS, polyLS, log4LS_acc)

            fit_rows.append({
                "run": str(run_root), "source": source_label, "mode": mode_suffix(normalize),
                "space": "LOGLOG", "x": "logS", "y": "logSig",
                "lin_slope": linLS["slope"], "lin_intercept": linLS["intercept"], "lin_r2": linLS["r2"], "lin_p": linLS["p"], "lin_aic": linLS["aic"],
                "poly_deg": poly_deg if poly_deg >= 2 else 0, "poly_r2": polyLS["r2"],
                "log4_attempted": bool(do_logistic4), "log4_accepted": bool(log4LS_acc is not None),
                "log4_status": log4LS_status, "log4_reason": log4LS_reason,
                "log4_L": (log4LS_acc["L"] if log4LS_acc else np.nan),
                "log4_k": (log4LS_acc["k"] if log4LS_acc else np.nan),
                "log4_x0": (log4LS_acc["x0"] if log4LS_acc else np.nan),
                "log4_c": (log4LS_acc["c"] if log4LS_acc else np.nan),
                "log4_r2": (log4LS_acc["r2"] if log4LS_acc else np.nan),
                "log4_p": (log4LS_acc["p_model"] if log4LS_acc else np.nan),
                "log4_aic": (log4LS_acc["aic"] if log4LS_acc else np.nan),
            })

        # Write NEW joined CSV per mode (so UNORM doesn't overwrite legacy)
        csv_new = results_dir / (f"joined_{tag}Hz.csv" if normalize else f"joined__UNORM_{tag}Hz.csv")
        df.to_csv(csv_new, index=False)

        # NEW plots (with suffix naming; legacy plots already created above for NORM)
        # Conductivity vs S (REAL and LOGLOG)
        out_real_sig = fname("cond_vs_S", normalize, False, tag, extra="LINEAR+POLY+LOG4" if do_logistic4 else "LINEAR+POLY")
        write_plot_xy(results_dir, csv_new.name, out_real_sig, source_label,
                      title=f"Conductivity vs Saturation (REAL) [{mode_suffix(normalize)}]",
                      xlabel="S", ylabel="sigma_imag",
                      x_col="S_real", y_col="sig_real",
                      y_lin_col="sig_lin_real",
                      y_poly_col=("sig_poly_real" if poly_deg >= 2 else None),
                      y_log4_col=("sig_log4_real" if do_logistic4 else None))

        out_log_sig = fname("cond_vs_S", normalize, True, tag, extra="LINEAR+POLY+LOG4" if do_logistic4 else "LINEAR+POLY")
        write_plot_xy(results_dir, csv_new.name, out_log_sig, source_label,
                      title=f"log10(sigma_imag) vs log10(S) [{mode_suffix(normalize)}]",
                      xlabel="logS", ylabel="logSig",
                      x_col="logS", y_col="logSig",
                      y_lin_col="logSig_lin",
                      y_poly_col=("logSig_poly" if poly_deg >= 2 else None),
                      y_log4_col=("logSig_log4" if do_logistic4 else None))

        # Also rho vs S
        out_real_rho = fname("rho_vs_S", normalize, False, tag, extra="LINEAR+POLY+LOG4" if do_logistic4 else "LINEAR+POLY")
        write_plot_xy(results_dir, csv_new.name, out_real_rho, source_label,
                      title=f"Resistivity vs Saturation (REAL) [{mode_suffix(normalize)}]",
                      xlabel="S", ylabel="rho",
                      x_col="S_real", y_col="rho_real",
                      y_lin_col="rho_lin_real",
                      y_poly_col=("rho_poly_real" if poly_deg >= 2 else None),
                      y_log4_col=("rho_log4_real" if do_logistic4 else None))

        out_log_rho = fname("rho_vs_S", normalize, True, tag, extra="LINEAR+POLY+LOG4" if do_logistic4 else "LINEAR+POLY")
        write_plot_xy(results_dir, csv_new.name, out_log_rho, source_label,
                      title=f"log10(rho) vs log10(S) [{mode_suffix(normalize)}]",
                      xlabel="logS", ylabel="logRho",
                      x_col="logS", y_col="logRho",
                      y_lin_col="logRho_lin",
                      y_poly_col=("logRho_poly" if (poly_raw and poly_deg >= 2) else None),
                      y_log4_col=("logRho_log4" if do_logistic4 else None))

        if not dry_run:
            # run new .gp scripts written by write_plot_xy
            for outpng in [out_real_sig, out_log_sig, out_real_rho, out_log_rho]:
                gp = Path(outpng).stem + ".gp"
                _run_gnuplot(gp, cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)

        # Excel sheets per mode & space
        excel_sheets[f"data__{mode_suffix(normalize)}__ALL"] = df.copy()

        out_rows.append({
            "run": str(run_root),
            "source": source_label,
            "mode": mode_suffix(normalize),
            "status": "OK",
            "points": int(len(df)),
            "logbook_sheet": used_sheet,
            "logbook_overlap": overlap_cnt,
            "logbook_block_rows": block_rows,
            "per_run_excel": str(excel_path),
            "joined_csv": str(csv_new),
        })

    # Write per-run Excel once (contains both modes)
    with pd.ExcelWriter(excel_path, engine="openpyxl") as xw:
        for sh, dfx in excel_sheets.items():
            # Excel sheet name limit is 31
            xw_name = sh[:31]
            dfx.to_excel(xw, sheet_name=xw_name, index=False)
        if fit_rows:
            pd.DataFrame(fit_rows).to_excel(xw, sheet_name="fits_summary", index=False)

    return out_rows


# =============================================================================
# MAIN CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("top", help="Top folder containing run subfolders")
    ap.add_argument("--freq", type=float, default=0.01, help="Target frequency in Hz")
    ap.add_argument("--sample", type=str, default=None,
                    help="Optional substring filter on SIP sheet names (debug)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Write scripts/CSVs but do not run gnuplot")
    ap.add_argument("--min-points", type=int, default=4,
                    help="Minimum matched points to accept a run")
    ap.add_argument("--max-df", type=float, default=0.0,
                    help="Max allowed |f_found - f_target| (0 disables)")
    ap.add_argument("--no-logistic", action="store_true",
                    help="Disable LEGACY logistic fit/plots (raw logRho) in legacy pipeline")

    ap.add_argument("--logistic-min-points", type=int, default=6,
                    help="Min points needed to attempt logistic fits")
    ap.add_argument("--logistic-fast", action="store_true",
                    help="Use SciPy-free fast logistic search")
    ap.add_argument("--logistic-max-seconds", type=float, default=8.0,
                    help="Time cap for logistic fit per run")
    ap.add_argument("--gnuplot-timeout", type=int, default=30,
                    help="Timeout seconds per gnuplot call")
    ap.add_argument("--clean-results", action="store_true",
                    help="Delete files inside each run's Results/ before writing outputs.")

    ap.add_argument("--ici-mode", type=str, default="sigref_over_sig",
                    choices=["sigref_over_sig","sig_over_sigref"],
                    help=("GOAL/INDEX ICI definition (legacy). VN target is fixed to sigref/sig."))

    ap.add_argument("--poly-deg", type=int, default=0,
                    help="Polynomial degree for NEW fits/plots (>=2). 0 disables.")
    ap.add_argument("--poly-raw", action="store_true",
                    help="Also apply polynomial overlay to logRho vs logS (NEW loglog).")

    ap.add_argument("--both-modes", action="store_true",
                    help="Run BOTH NORM (legacy preserved) and UNORM (new suffix outputs)")

    ap.add_argument("--logistic4", action="store_true",
                    help="Enable NEW broad 4-parameter logistic fits (REAL and LOGLOG)")

    args = ap.parse_args()

    base = Path(args.top).expanduser().resolve()
    if not base.exists():
        print(f"Not found: {base}")
        raise SystemExit(2)

    roots = discover_run_roots_fast(base)
    print(f"Found {len(roots)} run folder(s).", flush=True)

    rows = []
    t0 = time.time()

    for i, run in enumerate(roots, 1):
        print(f"[{i}/{len(roots)}] START {run}", flush=True)
        ti = time.time()
        try:
            res_rows = process_run(
                run_root=run,
                target_freq=args.freq,
                sheet_filter=args.sample,
                dry_run=args.dry_run,
                min_points=args.min_points,
                max_df=args.max_df,
                do_logistic_legacy=(not args.no_logistic),
                logistic_min_points=args.logistic_min_points,
                logistic_fast=args.logistic_fast,
                logistic_max_seconds=args.logistic_max_seconds,
                gnuplot_timeout=args.gnuplot_timeout,
                clean_results=args.clean_results,
                ici_mode=args.ici_mode,
                both_modes=args.both_modes,
                poly_deg=args.poly_deg,
                poly_raw=args.poly_raw,
                do_logistic4=args.logistic4,
            )
            rows.extend(res_rows)
            dt = time.time() - ti
            ok_modes = ",".join([r.get("mode","?") for r in res_rows if r.get("status") == "OK"])
            pts = res_rows[0].get("points", 0) if res_rows else 0
            print(f"[{i}/{len(roots)}] OK modes=[{ok_modes}] points={pts} ({dt:.1f}s)", flush=True)
        except subprocess.TimeoutExpired:
            rows.append({"run": str(run), "status": "FAIL", "reason": "gnuplot timeout", "points": 0})
            print(f"[{i}/{len(roots)}] FAIL gnuplot timeout", flush=True)
        except Exception as e:
            rows.append({"run": str(run), "status": "FAIL", "reason": str(e), "points": 0})
            print(f"[{i}/{len(roots)}] FAIL {e}", flush=True)

    tag = tag_freq(args.freq)
    summary_csv = base / f"Results_Summary_{tag}Hz.csv"
    summary_xlsx = base / f"Results_Summary_{tag}Hz.xlsx"

    if rows:
        df_sum = pd.DataFrame(rows)
        df_sum.to_csv(summary_csv, index=False)
        with pd.ExcelWriter(summary_xlsx, engine="openpyxl") as xw:
            df_sum.to_excel(xw, sheet_name="summary", index=False)

        print(f"\nSummary written: {summary_csv}", flush=True)
        print(f"Summary written: {summary_xlsx}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)

if __name__ == "__main__":
    main()
