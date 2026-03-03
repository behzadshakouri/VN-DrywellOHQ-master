#!/usr/bin/env python3
"""
auto_sip_pipeline_v21_noscipy.py

v21 (FIX + target plots):
- Keep RAW linear fits + optional logistic + separate logistic-only + diagnostics.
- Keep INDEX "goal" plot using through-origin fits on:
    logRI  = log10(rho/rho_ref)    = logRho - logRho_ref
    logICI = log10(sig_ref/sigma)  = logSig_ref - logImag   (FIXED)

- FIX v20 "went wrong" issue:
  * DO NOT hardcode numeric columns like using 6:11.
  * Instead, read the CSV header in Python and compute 1-based column indices,
    then write gnuplot scripts using those indices. This stays gnuplot-safe
    while remaining robust to any column reordering.

- Add "target" plots at selected frequency (publication-style, like your goals):
  * Saturation Exponent (Resistivity) at f: logRI vs logS with dashed blue fit
  * Saturation Exponent (Imag Conductivity Index) at f: logICI vs logS dashed red fit

- Add function to delete files inside Results folder (optional flag).

No SciPy required.
"""

import math
import argparse
import subprocess
import time
from pathlib import Path
import shutil

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

# Logistic constraints (avoid flat fits / overflow)
K_MIN = 0.05
K_MAX = 50.0

# Logistic acceptance settings
LOGISTIC_R2_EPS = 0.005
LOGISTIC_REQUIRE_BETTER_AIC = True


# ---------------- Small utilities ----------------
def _norm(s: str) -> str:
    return str(s).strip().lower()

def find_first_matching_file(folder: Path, hint: str, ext=".xlsx"):
    if not folder.exists():
        return None
    for p in folder.rglob(f"*{ext}"):
        if hint.lower() in p.name.lower():
            return p
    return None

def pick_col(df: pd.DataFrame, cands):
    low = {_norm(c): c for c in df.columns}
    for k in cands:
        kk = _norm(k)
        if kk in low:
            return low[kk]
    return ""

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


# ---------------- Results folder cleanup ----------------
def clear_results_folder(results_dir: Path, remove_subdirs: bool = False):
    """
    Deletes files inside Results folder. By default keeps subdirectories (if any).
    Safe: does not delete the Results folder itself.

    remove_subdirs=False: delete only files directly under results_dir
    remove_subdirs=True : delete everything under results_dir (files + subdirs)
    """
    if not results_dir.exists():
        return

    if remove_subdirs:
        for p in results_dir.iterdir():
            try:
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
            except Exception:
                pass
        return

    for p in results_dir.iterdir():
        try:
            if p.is_file():
                p.unlink()
        except Exception:
            pass


# ---------------- F distribution SF (no SciPy) ----------------
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

def normal_cdf(z):
    return 0.5*(1.0 + math.erf(z/math.sqrt(2.0)))


# ---------------- AIC ----------------
def aic_gaussian(ss_res, n, k_params):
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)


# ---------------- Fits ----------------
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


# ---------------- Logistic ----------------
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
        return True, f"R2_improved({lin_r2:.3f}->{log_r2:.3f})"
    if LOGISTIC_REQUIRE_BETTER_AIC and np.isfinite(lin_aic) and (log_aic < lin_aic):
        return True, f"AIC_improved({lin_aic:.2f}->{log_aic:.2f})"
    return False, f"rejected(R2 {lin_r2:.3f}->{log_r2:.3f}, AIC {lin_aic:.2f}->{log_aic:.2f})"


# ---------------- Logbook reading ----------------
def read_logbook_anysheet(log_xlsx: Path, sample_filter: str|None):
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

    sat_map = {}
    for _, r in df.iterrows():
        mid = str(r.get(id_col, "")).strip()
        if not mid or mid.lower() in ["nan","none"]:
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
            sat_map[mid] = float(S)

    if not sat_map:
        raise RuntimeError(f"Selected logbook sheet '{sh}' but mapped 0 saturation rows.")

    return sat_map, sh


# ---------------- SIP output reading ----------------
def read_sip_output(sip_xlsx: Path, target_freq: float, max_df: float):
    xls = pd.ExcelFile(sip_xlsx)
    out = []
    for sh in xls.sheet_names:
        df = pd.read_excel(sip_xlsx, sheet_name=sh)
        if df.empty:
            continue

        freq_col = pick_col(df, SIP_FREQ_CANDS) or df.columns[0]
        rho_col  = pick_col(df, SIP_RHO_CANDS)
        imag_col = pick_col(df, SIP_IMAG_CANDS)

        # Fallback: common xlsx layout
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
            "sigma_imag": float(sim)
        })
    return pd.DataFrame(out)


# ---------------- CSV header -> gnuplot numeric columns ----------------
def csv_col_map_1based(csv_path: Path):
    """
    Returns dict {colname: 1-based index} for gnuplot.
    Robust to any column order as long as names exist.
    """
    first = csv_path.read_text().splitlines()[0].strip()
    cols = [c.strip() for c in first.split(",")]
    return {name: (i+1) for i, name in enumerate(cols)}


# ---------------- Gnuplot runner ----------------
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


# ---------------- Gnuplot script writers (v21: header-mapped numeric cols) ----------------
def write_plot_index(results_dir: Path, csv_name: str, col, out_png: str,
                     target_freq: float, a_ri, r2_ri, n_idx, a_ici, r2_ici, p_idx):
    label = f"{target_freq:g}"
    x  = col["logS"]
    y1 = col["logRhoIdx"]
    y2 = col["logImagIdx"]

    gp = f"""reset
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
a1 = {a_ri:.12g}
f1(x) = a1*x
plot '{csv_name}' using {x}:{y1} with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('{label}Hz: n=%.2f, R^2=%.3f', {n_idx:.12g}, {r2_ri:.12g})

set title 'Imaginary Conductivity Index vs Saturation (log-log) at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
a2 = {a_ici:.12g}
f2(x) = a2*x
plot '{csv_name}' using {x}:{y2} with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('{label}Hz: p=%.2f, R^2=%.3f', {p_idx:.12g}, {r2_ici:.12g})

unset multiplot
"""
    (results_dir / "plot_index_exponents.gp").write_text(gp)


def write_plot_raw(results_dir: Path, csv_name: str, col, out_png: str, target_freq: float,
                   lin_rho, lin_im, log_fit_acc):
    label = f"{target_freq:g}"
    a_r = lin_rho["slope"]; b_r = lin_rho["intercept"]; r2_r = lin_rho["r2"]; p_r = lin_rho["p"]
    a_i = lin_im["slope"];  b_i = lin_im["intercept"];  r2_i = lin_im["r2"];  p_i = lin_im["p"]

    x  = col["logS"]
    yr = col["logRho"]
    yi = col["logImag"]

    log_defs = ""
    log_line = ""
    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        r2_log = log_fit_acc.get("r2", np.nan)
        p_model = log_fit_acc.get("p_model", np.nan)
        log_defs = (
            f"L1={L:.12g}\n"
            f"k1={k:.12g}\n"
            f"x01={x0:.12g}\n"
            f"c1={c:.12g}\n"
            "fL(x)=c1+L1/(1+exp(-k1*(x-x01)))\n"
        )
        if np.isfinite(r2_log) and np.isfinite(p_model):
            log_line = f", \\\n     fL(x) with lines dt 1 lw 4 title sprintf('logistic: R^2=%.3f, p=%.3g', {r2_log:.12g}, {p_model:.12g})"
        else:
            log_line = ", \\\n     fL(x) with lines dt 1 lw 4 title 'logistic fit'"

    gp = f"""reset
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
plot '{csv_name}' using {x}:{yr} with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('linear: n=%.2f, R^2=%.3f, p=%.3g', -a1, {r2_r:.12g}, {p_r:.12g}){log_line}

set title 'Imaginary Conductivity vs Saturation (log-log) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Imag Conductivity)'
a2={a_i:.12g}
b2={b_i:.12g}
f2(x)=a2*x+b2
plot '{csv_name}' using {x}:{yi} with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('linear: p=%.2f, R^2=%.3f, p=%.3g', a2, {r2_i:.12g}, {p_i:.12g})

unset multiplot
"""
    (results_dir / "plot_raw_exponents.gp").write_text(gp)


def write_plot_logistic_only(results_dir: Path, csv_name: str, col, out_png: str, target_freq: float, log_fit_acc):
    label = f"{target_freq:g}"
    x  = col["logS"]
    yr = col["logRho"]

    log_defs = ""
    log_line = ""
    if log_fit_acc is not None:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        r2_log = log_fit_acc.get("r2", np.nan)
        p_model = log_fit_acc.get("p_model", np.nan)
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

    gp = f"""reset
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
plot '{csv_name}' using {x}:{yr} with points pt 7 ps 1.7 title 'data', \\
     {log_line}
"""
    (results_dir / "plot_logistic_only.gp").write_text(gp)


def write_plot_logistic_diag(results_dir: Path, csv_name: str, col, out_png: str, target_freq: float, lin_rho, log_fit_acc):
    label = f"{target_freq:g}"
    a = lin_rho["slope"]; b = lin_rho["intercept"]
    has_log = (log_fit_acc is not None)

    x  = col["logS"]
    yr = col["logRho"]

    log_defs = ""
    if has_log:
        L = log_fit_acc["L"]; k = log_fit_acc["k"]; x0 = log_fit_acc["x0"]; c = log_fit_acc["c"]
        log_defs = (
            f"L1={L:.12g}\n"
            f"k1={k:.12g}\n"
            f"x01={x0:.12g}\n"
            f"c1={c:.12g}\n"
            "fL(x)=c1+L1/(1+exp(-k1*(x-x01)))\n"
        )

    panel1 = f"""plot '{csv_name}' using {x}:{yr} with points pt 7 ps 1.6 title 'data', \\
     f1(x) with lines dt 2 lw 3 title 'linear'"""
    if has_log:
        panel1 += """, \\
     fL(x) with lines dt 1 lw 4 title 'logistic (accepted)'"""
    panel1 += "\n"

    # IMPORTANT: residual panels use numeric columns with $<col>
    panel2 = f"plot '{csv_name}' using {x}:((${yr}) - f1(${x})) with points pt 7 ps 1.4 title 'residuals'\n"

    if has_log:
        panel3 = f"plot '{csv_name}' using {x}:((${yr}) - fL(${x})) with points pt 7 ps 1.4 title 'residuals'\n"
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

    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,1800 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 3,1

a1={a:.12g}
b1={b:.12g}
f1(x)=a1*x+b1
{log_defs}

set title 'Diagnostics: logRho vs logS at {label} Hz'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
{panel1}

set title 'Residuals: linear (data - linear)'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'residual'
{panel2}

set title 'Residuals: logistic (data - logistic)'
set grid
set xlabel 'log10(Saturation)'
set ylabel 'residual'
{extra3}
{panel3}

unset multiplot
"""
    (results_dir / "plot_logistic_diagnostics.gp").write_text(gp)


# --- Target-style exponent plots (single panel, blue n / red p) ---
def write_plot_target_exponent(results_dir: Path, csv_name: str, col, out_png: str,
                               title: str, xlab: str, ylab: str,
                               xcol_name: str, ycol_name: str,
                               slope_a: float, exp_val: float, r2: float,
                               color_name: str):
    """
    Single-panel target-style plot with dashed fit and legend:
      "<freq>Hz: n=..., R^2=..."
    color_name: 'blue' or 'red' (gnuplot color names)
    """
    x = col[xcol_name]
    y = col[ycol_name]

    # Make it look like your goal figure: grid, thicker line, larger points.
    gp = f"""reset
set datafile separator ','
set datafile missing ''
set term pngcairo size 1200,900 enhanced font 'Arial,26'
set output '{out_png}'

set grid
set key top right
set title '{title}'
set xlabel '{xlab}'
set ylabel '{ylab}'

a1 = {slope_a:.12g}
f1(x) = a1*x

set style line 1 lc rgb '{color_name}' lt 1 lw 3 dt 2
set style line 2 lc rgb '{color_name}' pt 7 ps 1.8

plot '{csv_name}' using {x}:{y} with points ls 2 title 'data', \\
     f1(x) with lines ls 1 title sprintf('%.2gHz: exp=%.2f, R^2=%.3f', {float('nan')}, {exp_val:.12g}, {r2:.12g})
"""
    # Replace the %.2gHz placeholder safely by writing a constant inside sprintf is messy;
    # we’ll just inject a fixed label via set key title-style alternative:
    # Instead, we’ll use a direct title string without sprintf frequency:
    gp = gp.replace(
        "title sprintf('%.2gHz: exp=%.2f, R^2=%.3f', nan,",
        "title sprintf('exp=%.2f, R^2=%.3f',"
    )
    (results_dir / f"plot_target_{out_png.replace('.png','')}.gp").write_text(gp)


def discover_run_roots_fast(base: Path):
    roots = set()
    for a in base.rglob(ANALYSIS_DIR):
        if not a.is_dir():
            continue
        run = a.parent
        if (run / RAW_DIR).is_dir():
            roots.add(run)
    return sorted(roots)


def process_run(run_root: Path, target_freq: float, sample_filter: str|None,
                dry_run: bool, min_points: int, max_df: float,
                do_logistic: bool, logistic_min_points: int,
                logistic_fast: bool, logistic_max_seconds: float,
                gnuplot_timeout: int,
                clear_results: bool,
                clear_results_all: bool):
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return None

    results_dir.mkdir(parents=True, exist_ok=True)
    if clear_results:
        clear_results_folder(results_dir, remove_subdirs=clear_results_all)

    sat_map, used_sheet = read_logbook_anysheet(log_xlsx, sample_filter)
    sip_df = read_sip_output(sip_xlsx, target_freq, max_df=max_df)

    tag = str(target_freq).replace(".","p")

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

    # Index (FIXED)
    sip_df["rho_idx"]     = sip_df["rho"] / rho_ref
    sip_df["imag_idx"]    = sig_ref / sip_df["sigma_imag"]          # ICI = sig_ref / sig
    sip_df["logRhoIdx"]   = sip_df["logRho"] - logRho_ref
    sip_df["logImagIdx"]  = logSig_ref - sip_df["logImag"]          # logICI = log(sig_ref) - log(sig)

    sip_df = sip_df[np.isfinite(sip_df["logRhoIdx"]) & np.isfinite(sip_df["logImagIdx"])]

    csv_path = results_dir / f"joined_{tag}Hz.csv"
    sip_df[[
        "measurement","freq_hz","S","rho","sigma_imag",
        "logS","logRho","logImag",
        "rho_idx","imag_idx","logRhoIdx","logImagIdx"
    ]].to_csv(csv_path, index=False)

    # Column map for gnuplot (v21 fix)
    col = csv_col_map_1based(csv_path)
    required = ["logS","logRho","logImag","logRhoIdx","logImagIdx"]
    for k in required:
        if k not in col:
            raise RuntimeError(f"CSV missing column '{k}' in {csv_path}")

    # RAW linear fits
    slope_r, intercept_r, r2_r, p_r, ss_res_lin_r, aic_lin_r = linear_fit(sip_df["logS"], sip_df["logRho"])
    slope_i, intercept_i, r2_i, p_i, ss_res_lin_i, aic_lin_i = linear_fit(sip_df["logS"], sip_df["logImag"])

    lin_rho = {"slope": slope_r, "intercept": intercept_r, "r2": r2_r, "p": p_r, "ss_res": ss_res_lin_r, "aic": aic_lin_r}
    lin_im  = {"slope": slope_i, "intercept": intercept_i, "r2": r2_i, "p": p_i, "ss_res": ss_res_lin_i, "aic": aic_lin_i}

    # INDEX through-origin fits (goal plot)
    a_ri,  r2_ri,  ss_ri  = linear_fit_through_origin(sip_df["logS"], sip_df["logRhoIdx"])
    a_ici, r2_ici, ss_ici = linear_fit_through_origin(sip_df["logS"], sip_df["logImagIdx"])
    n_idx = -a_ri
    p_idx = -a_ici

    # Logistic fit on RAW logRho vs logS
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

    # Report
    report = results_dir / f"fit_report_{tag}Hz.txt"
    report.write_text(
        f"Run root: {run_root}\n"
        f"SIP file: {sip_xlsx}\n"
        f"Logbook : {log_xlsx}\n"
        f"Logbook sheet used: {used_sheet}\n"
        f"Target freq requested: {target_freq} Hz\n"
        f"Matched points used: {len(sip_df)}\n"
        f"Reference point (max S): measurement={sip_df.iloc[iref]['measurement']}  S_ref={sip_df.iloc[iref]['S']:.6g}\n"
        f"rho_ref={rho_ref:.6g}, sigma_ref={sig_ref:.6g}\n\n"
        f"--- INDEX definitions (FIXED) ---\n"
        f"RI  = rho / rho_ref\n"
        f"ICI = sigma_ref / sigma_imag\n"
        f"logRI  = logRho  - logRho_ref\n"
        f"logICI = logSig_ref - logImag\n\n"
        f"--- INDEX through-origin fits ---\n"
        f"a_ri = {a_ri:.6g}, n(index) = {n_idx:.6g}, R2={r2_ri:.6g}, SSE={ss_ri:.6g}\n"
        f"a_ici= {a_ici:.6g}, p(index) = {p_idx:.6g}, R2={r2_ici:.6g}, SSE={ss_ici:.6g}\n\n"
        f"--- RAW Linear fit (logRho vs logS): logRho = a*logS + b ---\n"
        f"a = {slope_r:.6g}\n"
        f"b = {intercept_r:.6g}\n"
        f"n(raw) = {-slope_r:.6g}\n"
        f"R2 = {r2_r:.6g}\n"
        f"p(model) = {p_r:.6g}\n"
        f"SSE = {ss_res_lin_r:.6g}\n"
        f"AIC = {aic_lin_r:.6g}\n\n"
        f"--- RAW Linear fit (logImag vs logS): logImag = p*logS + q ---\n"
        f"p(raw) = {slope_i:.6g}\n"
        f"q = {intercept_i:.6g}\n"
        f"R2 = {r2_i:.6g}\n"
        f"p(model) = {p_i:.6g}\n"
        f"SSE = {ss_res_lin_i:.6g}\n"
        f"AIC = {aic_lin_i:.6g}\n\n"
        f"Logistic: status={log_status}, accepted={bool(log_fit_acc is not None)}, reason={log_accept_reason}\n"
    )

    # Standard plots (like v19)
    out_png_index = f"index_exponents_{tag}Hz.png"
    write_plot_index(results_dir, csv_path.name, col, out_png_index, target_freq,
                     a_ri, r2_ri, n_idx, a_ici, r2_ici, p_idx)

    out_png_raw = f"exponents_raw_{tag}Hz.png"
    write_plot_raw(results_dir, csv_path.name, col, out_png_raw, target_freq,
                   lin_rho, lin_im, log_fit_acc)

    out_png_log = f"exponents_logistic_only_{tag}Hz.png"
    write_plot_logistic_only(results_dir, csv_path.name, col, out_png_log, target_freq, log_fit_acc)

    out_png_diag = f"exponents_logistic_diagnostics_{tag}Hz.png"
    write_plot_logistic_diag(results_dir, csv_path.name, col, out_png_diag, target_freq, lin_rho, log_fit_acc)

    # Target-style plots (your goal look)
    # (These are based on INDEX definitions, which is what you want for Archie-style exponents.)
    sample_name = run_root.name.replace("_", " ")
    title_r = f"{sample_name}: Saturation Exponent (Resistivity Index) at {target_freq:g} Hz"
    title_i = f"{sample_name}: Saturation Exponent (Imag Conductivity Index) at {target_freq:g} Hz"

    # Use the same through-origin model: y = a*x, exponent = -a
    # Note: For RI plot, y=logRhoIdx; for ICI plot, y=logImagIdx.
    write_plot_target_exponent(
        results_dir, csv_path.name, col,
        out_png=f"TARGET_RI_exponent_{tag}Hz.png",
        title=title_r,
        xlab="log10(Saturation)",
        ylab="log10(Resistivity Index)",
        xcol_name="logS",
        ycol_name="logRhoIdx",
        slope_a=a_ri,
        exp_val=n_idx,
        r2=r2_ri,
        color_name="blue"
    )
    write_plot_target_exponent(
        results_dir, csv_path.name, col,
        out_png=f"TARGET_ICI_exponent_{tag}Hz.png",
        title=title_i,
        xlab="log10(Saturation)",
        ylab="log10(Imag Conductivity Index)",
        xcol_name="logS",
        ycol_name="logImagIdx",
        slope_a=a_ici,
        exp_val=p_idx,
        r2=r2_ici,
        color_name="red"
    )

    # Run gnuplot
    if not dry_run:
        _run_gnuplot("plot_index_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_raw_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_logistic_only.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_logistic_diagnostics.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

        # Target scripts have dynamic names
        _run_gnuplot(f"plot_target_TARGET_RI_exponent_{tag}Hz.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)
        _run_gnuplot(f"plot_target_TARGET_ICI_exponent_{tag}Hz.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

    return {
        "run": str(run_root),
        "status": "OK",
        "reason": "",
        "points": int(len(sip_df)),
        "logbook_sheet": used_sheet,
        "logistic_attempted": bool(log_fit is not None),
        "logistic_accepted": bool(log_fit_acc is not None),
        "logistic_status": log_status,
        "logistic_reason": log_accept_reason,
        "n_index": float(n_idx) if np.isfinite(n_idx) else np.nan,
        "p_index": float(p_idx) if np.isfinite(p_idx) else np.nan,
        "R2_RI": float(r2_ri) if np.isfinite(r2_ri) else np.nan,
        "R2_ICI": float(r2_ici) if np.isfinite(r2_ici) else np.nan,
    }


def main():
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

    # NEW: cleanup flags
    ap.add_argument("--clear-results", action="store_true",
                    help="Delete files inside each Results folder before writing outputs.")
    ap.add_argument("--clear-results-all", action="store_true",
                    help="If set with --clear-results, also delete subfolders inside Results.")

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
            res = process_run(
                run,
                target_freq=args.freq,
                sample_filter=args.sample,
                dry_run=args.dry_run,
                min_points=args.min_points,
                max_df=args.max_df,
                do_logistic=(not args.no_logistic),
                logistic_min_points=args.logistic_min_points,
                logistic_fast=args.logistic_fast,
                logistic_max_seconds=args.logistic_max_seconds,
                gnuplot_timeout=args.gnuplot_timeout,
                clear_results=args.clear_results,
                clear_results_all=args.clear_results_all
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

    tag = str(args.freq).replace(".","p")
    summary_path = base / f"Results_Summary_{tag}Hz.csv"
    if rows:
        pd.DataFrame(rows).to_csv(summary_path, index=False)
        print(f"\nSummary written: {summary_path}", flush=True)

    print(f"Done in {(time.time()-t0):.1f}s.", flush=True)


if __name__ == "__main__":
    main()
