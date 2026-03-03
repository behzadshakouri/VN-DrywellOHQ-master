#!/usr/bin/env python3
"""
auto_sip_pipeline_v25_noscipy.py

v25:
- Fix v24 SyntaxError: no f-string expression contains backslashes.
- Python 3.8+ compatible typing (no `str|None`).
- --clean-results removes files AND subfolders inside Results/ (keeps Results/ itself).
- TARGET plot: ONE figure (2 panels) using INDEX/TARGET definitions:
    RI  = rho/rho_ref,  logRI  = logRho - logRho_ref
    ICI = sigma_ref/sigma_imag, logICI = logSig_ref - logSig_imag   (FIXED)
  Through-origin "goal" fits.
- RAW plots unchanged (intercepted fits), logistic optional.
- Add RAWLIKE 2-panel figure (RAW resistivity + RAW imag) in same stacked style
  to match your "Van Nuys S12 Saturation Exponents..." example.
"""

import math
import argparse
import subprocess
import time
import shutil
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

import numpy as np
import pandas as pd

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


def _norm(s: str) -> str:
    return str(s).strip().lower()

def find_first_matching_file(folder: Path, hint: str, ext: str = ".xlsx") -> Optional[Path]:
    if not folder.exists():
        return None
    for p in folder.rglob(f"*{ext}"):
        if hint.lower() in p.name.lower():
            return p
    return None

def pick_col(df: pd.DataFrame, cands: List[str]) -> str:
    low = {_norm(c): c for c in df.columns}
    for k in cands:
        kk = _norm(k)
        if kk in low:
            return low[kk]
    return ""

def to_float(x: Any) -> float:
    try:
        if isinstance(x, str):
            xs = x.strip()
            if xs.endswith("%"):
                return float(xs[:-1].strip())
            return float(xs)
        return float(x)
    except Exception:
        return float("nan")

def safe_log10(x: Any) -> float:
    try:
        xv = float(x)
    except Exception:
        return float("nan")
    if not np.isfinite(xv) or xv <= 0:
        return float("nan")
    return math.log10(xv)

def sse(y: np.ndarray, yhat: np.ndarray) -> float:
    r = y - yhat
    return float(np.sum(r * r))

# ---------------- F distribution SF (no SciPy) ----------------
def betacf(a: float, b: float, x: float, maxit: int = 200, eps: float = 3e-14) -> float:
    am = 1.0; bm = 1.0; az = 1.0
    qab = a + b; qap = a + 1.0; qam = a - 1.0
    bz = 1.0 - qab * x / qap
    for m in range(1, maxit + 1):
        em = float(m); tem = em + em
        d = em * (b - em) * x / ((qam + tem) * (a + tem))
        ap = az + d * am
        bp = bz + d * bm
        d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem))
        app = ap + d * az
        bpp = bp + d * bz
        aold = az
        am = ap / bpp
        bm = bp / bpp
        az = app / bpp
        bz = 1.0
        if abs(az - aold) < eps * abs(az):
            return az
    return az

def betai(a: float, b: float, x: float) -> float:
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    ln_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    bt = math.exp(math.log(x) * a + math.log(1.0 - x) * b - ln_beta)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * betacf(a, b, x) / a
    else:
        return 1.0 - bt * betacf(b, a, 1.0 - x) / b

def f_dist_sf(F: float, d1: int, d2: int) -> float:
    if not np.isfinite(F) or F < 0:
        return float("nan")
    x = (d1 * F) / (d1 * F + d2)
    cdf = betai(d1 / 2.0, d2 / 2.0, x)
    return max(0.0, 1.0 - cdf)

def f_test_pvalue(ss_res: float, ss_tot: float, df_model: int, n: int) -> float:
    if not (np.isfinite(ss_res) and np.isfinite(ss_tot)):
        return float("nan")
    if ss_tot <= 0:
        return float("nan")
    if ss_res >= ss_tot:
        return 1.0
    df2 = n - df_model - 1
    if df2 <= 0:
        return float("nan")
    num = (ss_tot - ss_res) / df_model
    den = ss_res / df2
    if den <= 0:
        return float("nan")
    F = num / den
    return f_dist_sf(F, df_model, df2)

def aic_gaussian(ss_res: float, n: int, k_params: int) -> float:
    if n <= 0 or k_params <= 0 or (not np.isfinite(ss_res)) or ss_res <= 0:
        return float("nan")
    return float(n * math.log(ss_res / n) + 2.0 * k_params)

def linear_fit(x: Any, y: Any) -> Tuple[float,float,float,float,float,float]:
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), float("nan")
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    yhat = slope * x + intercept
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    pval = f_test_pvalue(ss_res, ss_tot, df_model=1, n=len(x))
    aic = aic_gaussian(ss_res, len(x), k_params=2)
    return float(slope), float(intercept), float(r2), float(pval), float(ss_res), float(aic)

def linear_fit_through_origin(x: Any, y: Any) -> Tuple[float,float,float]:
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan")
    den = float(np.dot(x, x))
    if den <= 0:
        return float("nan"), float("nan"), float("nan")
    a = float(np.dot(x, y) / den)
    yhat = a * x
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return a, float(r2), float(ss_res)

# ---------------- Logistic fit (fast grid search, no SciPy) ----------------
def logistic4(x: np.ndarray, L: float, k: float, x0: float, c: float) -> np.ndarray:
    z = k * (x - x0)
    z = np.clip(z, -60.0, 60.0)
    return c + L / (1.0 + np.exp(-z))

def _logistic_init_guess(x: np.ndarray, y: np.ndarray) -> Tuple[float,float,float,float]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    idx = np.argsort(x)
    x = x[idx]; y = y[idx]

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

def logistic_fit_noscipy(x: Any, y: Any, fast: bool = False, max_seconds: float = 0.0) -> Tuple[Optional[Dict[str,Any]], str]:
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

    stepL  = 0.25 * max(abs(Lb), 1e-6)
    stepk  = 0.25 * max(abs(kb), 1e-6)
    stepx0 = 0.20 * spanx
    stepc  = 0.25 * max(abs(Lb), 1e-6)

    def clamp_params(L: float, k: float, x0: float, c: float) -> Tuple[float,float,float,float]:
        L = max(1e-9, L)
        k = float(np.clip(k, K_MIN, K_MAX))
        return L, k, x0, c

    def current_sse(L: float, k: float, x0: float, c: float) -> float:
        return sse(y, logistic4(x, L, k, x0, c))

    Lb, kb, x0b, cb = clamp_params(Lb, kb, x0b, cb)
    best = current_sse(Lb, kb, x0b, cb)

    for _ in range(10):
        if max_seconds and (time.time() - t0) > max_seconds:
            return None, "timeout"
        improved = False
        for dL in [0, -stepL, stepL]:
            for dk in [0, -stepk, stepk]:
                for dx0 in [0, -stepx0, stepx0]:
                    for dc in [0, -stepc, stepc]:
                        if dL == dk == dx0 == dc == 0:
                            continue
                        if max_seconds and (time.time() - t0) > max_seconds:
                            return None, "timeout"
                        L2, k2, x02, c2 = clamp_params(Lb + dL, kb + dk, x0b + dx0, cb + dc)
                        s2 = current_sse(L2, k2, x02, c2)
                        if s2 + 1e-12 < best:
                            best = s2
                            Lb, kb, x0b, cb = L2, k2, x02, c2
                            improved = True
        stepL *= 0.6; stepk *= 0.6; stepx0 *= 0.6; stepc *= 0.6
        if (not improved) and max(stepL, stepk, stepx0, stepc) < 1e-4:
            break

    yhat = logistic4(x, Lb, kb, x0b, cb)
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    p_model = f_test_pvalue(ss_res, ss_tot, df_model=3, n=n)
    aic = aic_gaussian(ss_res, n, k_params=4)

    return {
        "L": float(Lb), "k": float(kb), "x0": float(x0b), "c": float(cb),
        "r2": float(r2),
        "p_model": float(p_model) if np.isfinite(p_model) else p_model,
        "aic": float(aic) if np.isfinite(aic) else aic,
        "ss_res": float(ss_res), "ss_tot": float(ss_tot),
    }, "ok"

def accept_logistic(lin_r2: float, lin_aic: float, log_fit: Optional[Dict[str,Any]]) -> Tuple[bool,str]:
    if log_fit is None:
        return False, "no_fit"
    log_r2 = log_fit.get("r2", float("nan"))
    log_aic = log_fit.get("aic", float("nan"))
    if (not np.isfinite(log_r2)) or (not np.isfinite(log_aic)):
        return False, "bad_metrics"
    if np.isfinite(lin_r2) and (log_r2 >= lin_r2 + LOGISTIC_R2_EPS):
        return True, "R2_improved"
    if LOGISTIC_REQUIRE_BETTER_AIC and np.isfinite(lin_aic) and (log_aic < lin_aic):
        return True, "AIC_improved"
    return False, "rejected"

# ---------------- Logbook reading ----------------
def read_logbook_anysheet(log_xlsx: Path, sample_filter: Optional[str]) -> Tuple[Dict[str,float], str]:
    xls = pd.ExcelFile(log_xlsx)
    best = None

    for sh in xls.sheet_names:
        try:
            df = pd.read_excel(log_xlsx, sheet_name=sh)
        except Exception:
            continue
        if df.empty or len(df.columns) < 2:
            continue

        id_col     = pick_col(df, LOGBOOK_ID_COLS_CAND)
        sat_col    = pick_col(df, LOGBOOK_SAT_COLS_CAND)
        sample_col = pick_col(df, LOGBOOK_SAMPLE_COLS_CAND)

        mw_col  = pick_col(df, LOGBOOK_MASSWATER_CAND)
        por_col = pick_col(df, LOGBOOK_POROSITY_CAND)
        vol_col = pick_col(df, LOGBOOK_COLVOL_CAND)

        score = 0
        if id_col: score += 3
        if sat_col: score += 4
        if sample_col: score += 1
        if (not sat_col) and (mw_col and por_col and vol_col): score += 3
        if score == 0:
            continue

        if (best is None) or (score > best[0]):
            best = (score, sh, df, id_col, sat_col, sample_col, mw_col, por_col, vol_col)

    if best is None:
        raise RuntimeError(f"Logbook has no sheet with saturation info. Sheets: {xls.sheet_names}")

    _, sh, df, id_col, sat_col, sample_col, mw_col, por_col, vol_col = best

    if sample_filter and sample_col:
        df = df[df[sample_col].astype(str).str.contains(sample_filter, case=False, na=False)]

    if not id_col:
        id_col = df.columns[0]

    sat_map: Dict[str,float] = {}
    for _, r in df.iterrows():
        mid = str(r.get(id_col, "")).strip()
        if not mid or mid.lower() in ["nan", "none"]:
            continue

        S = float("nan")
        if sat_col:
            S = to_float(r.get(sat_col, float("nan")))
            if np.isfinite(S) and S > 1.0 and S <= 100.0:
                S = S / 100.0
        else:
            Mw  = to_float(r.get(mw_col, float("nan")))   # g
            phi = to_float(r.get(por_col, float("nan")))
            V   = to_float(r.get(vol_col, float("nan")))  # m3
            if np.isfinite(Mw) and np.isfinite(phi) and np.isfinite(V) and phi > 0 and V > 0:
                Mw_kg = Mw / 1000.0
                rho_w = 1000.0
                S = Mw_kg / (rho_w * (phi * V))

        if np.isfinite(S) and S > 0:
            sat_map[mid] = float(S)

    if not sat_map:
        raise RuntimeError(f"Selected logbook sheet '{sh}' but mapped 0 saturation rows.")

    return sat_map, sh

# ---------------- SIP output reading ----------------
def read_sip_output(sip_xlsx: Path, target_freq: float, max_df: float) -> pd.DataFrame:
    xls = pd.ExcelFile(sip_xlsx)
    out: List[Dict[str,Any]] = []
    for sh in xls.sheet_names:
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
            "freq_hz": f_found,
            "rho": float(rho),
            "sigma_imag": float(sim),
        })
    return pd.DataFrame(out)

# ---------------- Results cleaning ----------------
def clean_results_dir(results_dir: Path) -> None:
    """Delete EVERYTHING inside Results/ (files + subfolders), keep Results/ itself."""
    if not results_dir.exists():
        return
    if results_dir.name != RESULTS_DIR:
        raise RuntimeError(f"Refusing to clean non-Results folder: {results_dir}")
    for p in results_dir.iterdir():
        try:
            if p.is_dir():
                shutil.rmtree(p)
            else:
                p.unlink()
        except Exception:
            pass

# ---------------- Plot writers (gnuplot; numeric cols only) ----------------
def write_plot_target(results_dir: Path, csv_name: str, out_png: str, run_name: str,
                      target_freq: float, n_idx: float, r2_ri: float, p_idx: float, r2_ici: float) -> None:
    label = f"{target_freq:g}"

    # Target uses INDEX columns:
    # x=logS (col6)
    # y=logRhoIdx (col11)  -> RI target
    # y=logImagIdx (col12) -> ICI target (FIXED)
    gp = (
f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set tics out

# --- TOP: RI (TARGET) ---
set title '{run_name}: Saturation Exponent (Resistivity Index) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity Index)'
set key top right

a1 = {-n_idx:.12g}
f1(x) = a1*x
plot '{csv_name}' using 6:11 with points pt 7 ps 1.8 lc rgb '#0066ff' title 'data', \\
     f1(x) with lines dt 2 lw 3 lc rgb '#0066ff' title sprintf('{label}Hz: n = %.2f, R^2 = %.3f', {n_idx:.12g}, {r2_ri:.12g})

# --- BOTTOM: ICI (TARGET) ---
set title '{run_name}: Saturation Exponent (Imag Conductivity Index) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Imag Conductivity Index)'
set key top right

a2 = {-p_idx:.12g}
f2(x) = a2*x
plot '{csv_name}' using 6:12 with points pt 7 ps 1.8 lc rgb '#dd0000' title 'data', \\
     f2(x) with lines dt 2 lw 3 lc rgb '#dd0000' title sprintf('{label}Hz: p = %.2f, R^2 = %.3f', {p_idx:.12g}, {r2_ici:.12g})

unset multiplot
"""
    )
    (results_dir / "plot_target_exponents.gp").write_text(gp)

def write_plot_rawlike_two_panel(results_dir: Path, csv_name: str, out_png: str, run_name: str,
                                 target_freq: float,
                                 slope_r: float, intercept_r: float, r2_r: float,
                                 slope_i: float, intercept_i: float, r2_i: float) -> None:
    """
    This is the "Van Nuys S12 Saturation Exponents ..." style you showed:
    - TOP: RAW logRho vs logS (intercepted)
    - BOTTOM: RAW logImag vs logS (intercepted)
    It is NOT TARGET.
    """
    label = f"{target_freq:g}"
    gp = (
f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set tics out

# --- TOP: RAW RHO ---
set title '{run_name}: Saturation Exponent (Resistivity) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
set key top right

a1={slope_r:.12g}
b1={intercept_r:.12g}
f1(x)=a1*x+b1
plot '{csv_name}' using 6:7 with points pt 7 ps 1.8 lc rgb '#0066ff' title 'data', \\
     f1(x) with lines dt 2 lw 3 lc rgb '#0066ff' title sprintf('{label}Hz: n = %.2f, R^2 = %.3f', -a1, {r2_r:.12g})

# --- BOTTOM: RAW IMAG ---
set title '{run_name}: Saturation Exponent (Imaginary Conductivity) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Imaginary Conductivity)'
set key top right

a2={slope_i:.12g}
b2={intercept_i:.12g}
f2(x)=a2*x+b2
plot '{csv_name}' using 6:8 with points pt 7 ps 1.8 lc rgb '#dd0000' title 'data', \\
     f2(x) with lines dt 2 lw 3 lc rgb '#dd0000' title sprintf('{label}Hz: p = %.2f, R^2 = %.3f', a2, {r2_i:.12g})

unset multiplot
"""
    )
    (results_dir / "plot_rawlike_exponents.gp").write_text(gp)

def write_plot_index(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                     n_idx: float, r2_ri: float, p_idx: float, r2_ici: float) -> None:
    label = f"{target_freq:g}"
    gp = (
f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set key top right
set xlabel 'log10(Saturation)'

set title 'Resistivity Index vs Saturation (log-log) at {label} Hz'
set ylabel 'log10(Resistivity Index)'
a1 = {-n_idx:.12g}
f1(x) = a1*x
plot '{csv_name}' using 6:11 with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('{label}Hz: n=%.2f, R^2=%.3f', {n_idx:.12g}, {r2_ri:.12g})

set title 'Imag Conductivity Index vs Saturation (log-log) at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
a2 = {-p_idx:.12g}
f2(x) = a2*x
plot '{csv_name}' using 6:12 with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('{label}Hz: p=%.2f, R^2=%.3f', {p_idx:.12g}, {r2_ici:.12g})

unset multiplot
"""
    )
    (results_dir / "plot_index_exponents.gp").write_text(gp)

def write_plot_raw(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
                   lin_rho: Dict[str,float], lin_im: Dict[str,float], log_fit_acc: Optional[Dict[str,Any]]) -> None:
    label = f"{target_freq:g}"
    a_r = lin_rho["slope"]; b_r = lin_rho["intercept"]; r2_r = lin_rho["r2"]; p_r = lin_rho["p"]; aic_r = lin_rho["aic"]
    a_i = lin_im["slope"];  b_i = lin_im["intercept"];  r2_i = lin_im["r2"];  p_i = lin_im["p"]

    log_defs = ""
    log_line = ""
    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        r2_log = log_fit_acc.get("r2", float("nan"))
        p_model = log_fit_acc.get("p_model", float("nan"))
        log_defs = (
            f"L1={L:.12g}\n"
            f"k1={k:.12g}\n"
            f"x01={x0:.12g}\n"
            f"c1={c:.12g}\n"
            "fL(x)=c1+L1/(1+exp(-k1*(x-x01)))\n"
        )
        if np.isfinite(r2_log) and np.isfinite(p_model):
            log_line = (
                ", \\\n"
                f"     fL(x) with lines dt 1 lw 4 title sprintf('logistic: R^2=%.3f, p=%.3g', {r2_log:.12g}, {p_model:.12g})"
            )
        else:
            log_line = ", \\\n     fL(x) with lines dt 1 lw 4 title 'logistic fit'"

    gp = (
f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'
set multiplot layout 2,1

set grid
set key top right

set title 'Resistivity vs Saturation (log-log) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
a1={a_r:.12g}
b1={b_r:.12g}
f1(x)=a1*x+b1
{log_defs}
plot '{csv_name}' using 6:7 with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('linear: n=%.2f, R^2=%.3f, p=%.3g', -a1, {r2_r:.12g}, {p_r:.12g}){log_line}

set title 'Imag Conductivity vs Saturation (log-log) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Imag Conductivity)'
a2={a_i:.12g}
b2={b_i:.12g}
f2(x)=a2*x+b2
plot '{csv_name}' using 6:8 with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('linear: p=%.2f, R^2=%.3f, p=%.3g', a2, {r2_i:.12g}, {p_i:.12g})

unset multiplot
"""
    )
    (results_dir / "plot_raw_exponents.gp").write_text(gp)

def write_plot_logistic_only(results_dir: Path, csv_name: str, out_png: str,
                             target_freq: float, log_fit_acc: Optional[Dict[str,Any]]) -> None:
    label = f"{target_freq:g}"
    log_defs = ""
    log_line = ""
    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        r2_log = log_fit_acc.get("r2", float("nan"))
        p_model = log_fit_acc.get("p_model", float("nan"))
        log_defs = (
            f"L1={L:.12g}\n"
            f"k1={k:.12g}\n"
            f"x01={x0:.12g}\n"
            f"c1={c:.12g}\n"
            "fL(x)=c1+L1/(1+exp(-k1*(x-x01)))\n"
        )
        if np.isfinite(r2_log) and np.isfinite(p_model):
            log_line = f"fL(x) with lines dt 1 lw 4 title sprintf('logistic: R^2=%.3f, p=%.3g', {r2_log:.12g}, {p_model:.12g})"
        else:
            log_line = "fL(x) with lines dt 1 lw 4 title 'logistic fit'"

    if not log_line:
        log_line = "1/0 notitle"

    gp = (
f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,800 enhanced font 'Arial,26'
set output '{out_png}'

set title 'Resistivity vs Saturation (log-log) at {label} Hz (logistic only)'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'

{log_defs}
plot '{csv_name}' using 6:7 with points pt 7 ps 1.7 title 'data', \\
     {log_line}
"""
    )
    (results_dir / "plot_logistic_only.gp").write_text(gp)

def write_plot_logistic_diag(results_dir: Path, csv_name: str, out_png: str,
                             target_freq: float, has_log: bool) -> None:
    label = f"{target_freq:g}"

    if has_log:
        # plot residuals in col14 (resid_log_rho)
        extra3 = ""
        panel3 = f"plot '{csv_name}' using 6:14 with points pt 7 ps 1.4 title 'residuals'\n"
        extra_log_curve = ", \\\n     '{0}' using 6:16 with lines dt 1 lw 4 title 'logistic (accepted)'".format(csv_name)
    else:
        extra_log_curve = ""
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

    gp = (
f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1800 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 3,1

set title 'Diagnostics: logRho vs logS at {label} Hz'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
plot '{csv_name}' using 6:7 with points pt 7 ps 1.6 title 'data', \\
     '{csv_name}' using 6:15 with lines dt 2 lw 3 title 'linear'{extra_log_curve}

set title 'Residuals: linear (data - linear)'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'residual'
plot '{csv_name}' using 6:13 with points pt 7 ps 1.4 title 'residuals'

set title 'Residuals: logistic (data - logistic)'
set grid
set xlabel 'log10(Saturation)'
set ylabel 'residual'
{extra3}
{panel3}

unset multiplot
"""
    )
    (results_dir / "plot_logistic_diagnostics.gp").write_text(gp)

def discover_run_roots_fast(base: Path) -> List[Path]:
    roots = set()
    for a in base.rglob(ANALYSIS_DIR):
        if not a.is_dir():
            continue
        run = a.parent
        if (run / RAW_DIR).is_dir():
            roots.add(run)
    return sorted(roots)

def _run_gnuplot(script_name: str, cwd: Path, timeout_s: int, fatal: bool) -> Tuple[bool,str]:
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

def process_run(run_root: Path, target_freq: float, sample_filter: Optional[str],
                dry_run: bool, min_points: int, max_df: float,
                do_logistic: bool, logistic_min_points: int,
                logistic_fast: bool, logistic_max_seconds: float,
                gnuplot_timeout: int,
                clean_results: bool) -> Optional[Dict[str,Any]]:

    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return None

    results_dir.mkdir(parents=True, exist_ok=True)
    if clean_results:
        clean_results_dir(results_dir)

    sat_map, used_sheet = read_logbook_anysheet(log_xlsx, sample_filter)
    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df)

    tag = str(target_freq).replace(".", "p")

    if sip_df.empty:
        return {"run": str(run_root), "status": "FAIL",
                "reason": "No SIP sheets readable / target freq missing",
                "points": 0, "logbook_sheet": used_sheet}

    sip_df["S"] = sip_df["measurement"].map(sat_map).astype(float)
    sip_df["logS"]    = sip_df["S"].apply(safe_log10)
    sip_df["logRho"]  = sip_df["rho"].apply(safe_log10)
    sip_df["logImag"] = sip_df["sigma_imag"].apply(safe_log10)

    sip_df = sip_df[np.isfinite(sip_df["logS"]) & np.isfinite(sip_df["logRho"]) & np.isfinite(sip_df["logImag"])]

    if len(sip_df) < min_points:
        return {"run": str(run_root), "status": "FAIL",
                "reason": f"Not enough matched points (need >= {min_points})",
                "points": int(len(sip_df)), "logbook_sheet": used_sheet}

    # Reference: max S
    iref = int(np.nanargmax(sip_df["S"].to_numpy(float)))
    rho_ref  = float(sip_df.iloc[iref]["rho"])
    sig_ref  = float(sip_df.iloc[iref]["sigma_imag"])
    logRho_ref = safe_log10(rho_ref)
    logSig_ref = safe_log10(sig_ref)

    # Index (FIXED ICI)
    sip_df["rho_idx"]     = sip_df["rho"] / rho_ref
    sip_df["imag_idx"]    = sig_ref / sip_df["sigma_imag"]
    sip_df["logRhoIdx"]   = sip_df["logRho"] - logRho_ref
    sip_df["logImagIdx"]  = logSig_ref - sip_df["logImag"]

    sip_df = sip_df[np.isfinite(sip_df["logRhoIdx"]) & np.isfinite(sip_df["logImagIdx"])]

    # RAW fits (intercepted)
    slope_r, intercept_r, r2_r, p_r, ss_res_lin_r, aic_lin_r = linear_fit(sip_df["logS"], sip_df["logRho"])
    slope_i, intercept_i, r2_i, p_i, ss_res_lin_i, aic_lin_i = linear_fit(sip_df["logS"], sip_df["logImag"])
    lin_rho = {"slope": slope_r, "intercept": intercept_r, "r2": r2_r, "p": p_r, "aic": aic_lin_r}
    lin_im  = {"slope": slope_i, "intercept": intercept_i, "r2": r2_i, "p": p_i, "aic": aic_lin_i}

    # TARGET/INDEX through-origin fits
    a_ri,  r2_ri,  _ = linear_fit_through_origin(sip_df["logS"], sip_df["logRhoIdx"])
    a_ici, r2_ici, _ = linear_fit_through_origin(sip_df["logS"], sip_df["logImagIdx"])
    n_idx = -a_ri
    p_idx = -a_ici

    # Logistic (on RAW logRho vs logS)
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

    # Precompute residual columns (gnuplot never needs expressions)
    logS = sip_df["logS"].to_numpy(float)
    logR = sip_df["logRho"].to_numpy(float)

    lin_pred = slope_r * logS + intercept_r
    sip_df["resid_lin_rho"] = logR - lin_pred

    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        log_pred = logistic4(logS, L, k, x0, c)
        sip_df["resid_log_rho"] = logR - log_pred
        sip_df["logRho_linfit"] = lin_pred
        sip_df["logRho_logfit"] = log_pred
    else:
        sip_df["resid_log_rho"] = np.nan
        sip_df["logRho_linfit"] = lin_pred
        sip_df["logRho_logfit"] = np.nan

    # Write CSV (fixed column order)
    csv_path = results_dir / f"joined_{tag}Hz.csv"
    sip_df[[
        "measurement","freq_hz","S","rho","sigma_imag",
        "logS","logRho","logImag",
        "rho_idx","imag_idx","logRhoIdx","logImagIdx",
        "resid_lin_rho","resid_log_rho",
        "logRho_linfit","logRho_logfit"
    ]].to_csv(csv_path, index=False)

    report = results_dir / f"fit_report_{tag}Hz.txt"
    report.write_text(
        "Run root: {run}\n"
        "Matched points used: {npts}\n"
        "Reference point (max S): {mref}  S_ref={sref:.6g}\n"
        "rho_ref={rref:.6g}, sigma_ref={srefsig:.6g}\n\n"
        "INDEX defs: RI=rho/rho_ref; ICI=sigma_ref/sigma_imag\n"
        "n(index)={nidx:.6g} (R2={r2ri:.6g})\n"
        "p(index)={pidx:.6g} (R2={r2ici:.6g})\n\n"
        "RAW linear logRho=a*logS+b: a={ar:.6g}, b={br:.6g}, n(raw)={nraw:.6g}, R2={r2r:.6g}\n"
        "RAW linear logImag=p*logS+q: p={ai:.6g}, q={bi:.6g}, R2={r2i:.6g}\n\n"
        "Logistic: status={ls}, accepted={acc}, reason={reason}\n"
        .format(
            run=run_root,
            npts=len(sip_df),
            mref=str(sip_df.iloc[iref]["measurement"]),
            sref=float(sip_df.iloc[iref]["S"]),
            rref=rho_ref,
            srefsig=sig_ref,
            nidx=n_idx, r2ri=r2_ri,
            pidx=p_idx, r2ici=r2_ici,
            ar=slope_r, br=intercept_r, nraw=-slope_r, r2r=r2_r,
            ai=slope_i, bi=intercept_i, r2i=r2_i,
            ls=log_status,
            acc=bool(log_fit_acc is not None),
            reason=log_accept_reason
        )
    )

    run_name = run_root.name

    # Plot scripts
    write_plot_target(results_dir, csv_path.name, f"target_exponents_{tag}Hz.png",
                      run_name, target_freq, n_idx, r2_ri, p_idx, r2_ici)

    write_plot_index(results_dir, csv_path.name, f"index_exponents_{tag}Hz.png",
                     target_freq, n_idx, r2_ri, p_idx, r2_ici)

    write_plot_raw(results_dir, csv_path.name, f"exponents_raw_{tag}Hz.png",
                   target_freq, lin_rho, lin_im, log_fit_acc)

    write_plot_rawlike_two_panel(results_dir, csv_path.name, f"rawlike_exponents_{tag}Hz.png",
                                 run_name, target_freq,
                                 slope_r, intercept_r, r2_r,
                                 slope_i, intercept_i, r2_i)

    write_plot_logistic_only(results_dir, csv_path.name, f"exponents_logistic_only_{tag}Hz.png",
                             target_freq, log_fit_acc)

    write_plot_logistic_diag(results_dir, csv_path.name, f"exponents_logistic_diagnostics_{tag}Hz.png",
                             target_freq, has_log=(log_fit_acc is not None))

    # Run gnuplot
    if not dry_run:
        _run_gnuplot("plot_target_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_index_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_raw_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_rawlike_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_logistic_only.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_logistic_diagnostics.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

    return {
        "run": str(run_root),
        "status": "OK",
        "points": int(len(sip_df)),
        "logistic_attempted": bool(log_fit is not None),
        "logistic_accepted": bool(log_fit_acc is not None),
        "logistic_status": log_status,
        "logistic_reason": log_accept_reason,
    }

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("top")
    ap.add_argument("--freq", type=float, default=0.01)
    ap.add_argument("--sample", type=str, default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--min-points", type=int, default=4)
    ap.add_argument("--max-df", type=float, default=0.0)
    ap.add_argument("--no-logistic", action="store_true")
    ap.add_argument("--logistic-min-points", type=int, default=6)
    ap.add_argument("--logistic-fast", action="store_true")
    ap.add_argument("--logistic-max-seconds", type=float, default=8.0)
    ap.add_argument("--gnuplot-timeout", type=int, default=30)
    ap.add_argument("--clean-results", action="store_true",
                    help="Delete files + subfolders inside each run's Results/ before writing new outputs.")
    args = ap.parse_args()

    base = Path(args.top).expanduser().resolve()
    if not base.exists():
        print(f"Not found: {base}")
        raise SystemExit(2)

    roots = discover_run_roots_fast(base)
    print(f"Found {len(roots)} run folder(s).", flush=True)

    rows: List[Dict[str,Any]] = []
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
                clean_results=args.clean_results
            )
            if res is None:
                continue
            rows.append(res)
            dt = time.time() - ti
            print(f"[{i}/{len(roots)}] OK points={res.get('points',0)} "
                  f"logistic_attempted={res.get('logistic_attempted',False)} "
                  f"logistic_accepted={res.get('logistic_accepted',False)} ({dt:.1f}s)", flush=True)
        except subprocess.TimeoutExpired:
            rows.append({"run": str(run), "status": "FAIL", "reason": "gnuplot timeout", "points": 0})
            print(f"[{i}/{len(roots)}] FAIL gnuplot timeout", flush=True)
        except Exception as e:
            rows.append({"run": str(run), "status": "FAIL", "reason": str(e), "points": 0})
            print(f"[{i}/{len(roots)}] FAIL {e}", flush=True)

    tag = str(args.freq).replace(".", "p")
    summary_path = base / f"Results_Summary_{tag}Hz.csv"
    if rows:
        pd.DataFrame(rows).to_csv(summary_path, index=False)
        print(f"\nSummary written: {summary_path}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)

if __name__ == "__main__":
    main()
