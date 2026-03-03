#!/usr/bin/env python3
"""
auto_sip_pipeline_v27_noscipy.py

v27 (FULL, corrected):
- TARGET plot now matches your "Van Nuys S12 Saturation Exponents at 0.01 Hz":
  * Uses RAW-like (WITH intercept) linear fits on VN-style index axes:
      logRI_vn  = log10(rho/rho_ref)        = logRho - logRho_ref
      logICI_vn = log10(sig_ref/sig_imag)   = logSig_ref - logImag   (always downtrend)
    -> TARGET legend shows n,p and R^2 that match RAW behavior (e.g., n≈1.03, p≈0.72)
    -> BOTH panels trend downward by construction.

- INDEX plot stays your "goal/through-origin" plot:
    logRhoIdx  = log10(rho/rho_ref)
    logImagIdx depends on --ici-mode:
      --ici-mode sigref_over_sig  => logImagIdx = log10(sig_ref/sig_imag)
      --ici-mode sig_over_sigref  => logImagIdx = log10(sig_imag/sig_ref)
  p(index) is reported as a POSITIVE exponent consistently.

- Uses robust gnuplot column("name") everywhere (no fragile column numbers).
- Keeps v26 features: no SciPy F-test/AIC, optional logistic on RAW logRho vs logS, diagnostics,
  --clean-results, etc.

Runner example:
  python3 auto_sip_pipeline_v27_noscipy.py "/mnt/3rd900/Projects/LA Project new/For_LBNL/SIP" \
    --freq 0.01 --min-points 4 \
    --logistic-min-points 6 --logistic-fast --logistic-max-seconds 8 \
    --clean-results --ici-mode sig_over_sigref
"""

import math
import argparse
import subprocess
import time
from pathlib import Path

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

def aic_gaussian(ss_res, n, k_params):
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)

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

# ---------------- Logistic fit (fast grid search, no SciPy) ----------------
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

# ---------------- Results cleaning ----------------
def clean_results_dir(results_dir: Path):
    """Delete files inside Results/ (keep folder; do not delete subdirs)."""
    if not results_dir.exists():
        return
    for p in results_dir.iterdir():
        try:
            if p.is_file() or p.is_symlink():
                p.unlink()
        except Exception:
            pass

# ---------------- Plot writers ----------------
def write_plot_target(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
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

# --- TOP: RI (Van Nuys raw-like, WITH intercept fit) ---
set title 'VNs-12: Saturation Exponent (Resistivity Index) at {label} Hz'
set ylabel 'log10(Resistivity Index)'
plot '{csv}' using (column("logS")):(column("logRI_vn")) with points pt 7 ps 1.8 lc rgb '#0066ff' title 'data', \\
     '{csv}' using (column("logS")):(column("logRI_vn_fit")) with lines dt 2 lw 3 lc rgb '#0066ff' \\
     title sprintf('{label}Hz: n = %.2f, R^2 = %.3f', {nval}, {r2a})

# --- BOTTOM: ICI (Van Nuys raw-like, WITH intercept fit; always downtrend) ---
set title 'VNs-12: Saturation Exponent (Imag Conductivity Index) at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
plot '{csv}' using (column("logS")):(column("logICI_vn")) with points pt 7 ps 1.8 lc rgb '#dd0000' title 'data', \\
     '{csv}' using (column("logS")):(column("logICI_vn_fit")) with lines dt 2 lw 3 lc rgb '#dd0000' \\
     title sprintf('{label}Hz: p = %.2f, R^2 = %.3f', {pval}, {r2b})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label,
           nval=n_vn, r2a=r2_ri_vn, pval=p_vn, r2b=r2_ici_vn)
    (results_dir / "plot_target_exponents.gp").write_text(gp)

def write_plot_index(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
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

# --- TOP: goal RI (through-origin) ---
set title 'Resistivity Index vs Saturation (log-log) at {label} Hz'
set ylabel 'log10(Resistivity Index)'
a1 = {a1}
f1(x) = a1*x
plot '{csv}' using (column("logS")):(column("logRhoIdx")) with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('{label}Hz: n=%.2f, R^2=%.3f', {n_idx}, {r2_ri})

# --- BOTTOM: goal ICI (through-origin) ---
set title 'Imag Conductivity Index vs Saturation (log-log) at {label} Hz'
set ylabel 'log10(Imag Conductivity Index)'
a2 = {a2}
f2(x) = a2*x
plot '{csv}' using (column("logS")):(column("logImagIdx")) with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('{label}Hz: p=%.2f, R^2=%.3f', {p_idx}, {r2_ici})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label,
           a1=(-n_idx), n_idx=n_idx, r2_ri=r2_ri,
           a2=(-p_idx), p_idx=p_idx, r2_ici=r2_ici)
    (results_dir / "plot_index_exponents.gp").write_text(gp)

def write_plot_raw(results_dir: Path, csv_name: str, out_png: str, target_freq: float,
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

set title 'Resistivity vs Saturation (log-log) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'
a1={a1:.12g}
b1={b1:.12g}
f1(x)=a1*x+b1
{log_defs}
plot '{csv}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('linear: n=%.2f, R^2=%.3f, p=%.3g', -a1, {r2_r:.12g}, {p_r:.12g}){log_clause}

set title 'Imag Conductivity vs Saturation (log-log) at {label} Hz'
set xlabel 'log10(Saturation)'
set ylabel 'log10(Imag Conductivity)'
a2={a2:.12g}
b2={b2:.12g}
f2(x)=a2*x+b2
plot '{csv}' using (column("logS")):(column("logImag")) with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('linear: p=%.2f, R^2=%.3f, p=%.3g', a2, {r2_i:.12g}, {p_i:.12g})

unset multiplot
""".format(out_png=out_png, csv=csv_name, label=label,
           a1=a_r, b1=b_r, r2_r=r2_r, p_r=p_r,
           a2=a_i, b2=b_i, r2_i=r2_i, p_i=p_i,
           log_defs=log_defs, log_clause=log_clause)
    (results_dir / "plot_raw_exponents.gp").write_text(gp)

def write_plot_logistic_only(results_dir: Path, csv_name: str, out_png: str, target_freq: float, log_fit_acc):
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

set title 'Resistivity vs Saturation (log-log) at {label} Hz (logistic only)'
set grid
set key top right
set xlabel 'log10(Saturation)'
set ylabel 'log10(Resistivity)'

{log_defs}
plot '{csv}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.7 title 'data', \\
     {log_line}
""".format(out_png=out_png, csv=csv_name, label=label, log_defs=log_defs, log_line=log_line)
    (results_dir / "plot_logistic_only.gp").write_text(gp)

def write_plot_logistic_diag(results_dir: Path, csv_name: str, out_png: str, target_freq: float, has_log: bool):
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

set title 'Diagnostics: logRho vs logS at {label} Hz'
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
""".format(out_png=out_png, label=label, csv=csv_name,
           extra_curve=extra_curve, extra3=extra3, panel3=panel3)
    (results_dir / "plot_logistic_diagnostics.gp").write_text(gp)

# ---------------- Runner helpers ----------------
def discover_run_roots_fast(base: Path):
    roots = set()
    for a in base.rglob(ANALYSIS_DIR):
        if not a.is_dir():
            continue
        run = a.parent
        if (run / RAW_DIR).is_dir():
            roots.add(run)
    return sorted(roots)

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

def process_run(run_root: Path, target_freq: float, sample_filter: str|None,
                dry_run: bool, min_points: int, max_df: float,
                do_logistic: bool, logistic_min_points: int,
                logistic_fast: bool, logistic_max_seconds: float,
                gnuplot_timeout: int,
                clean_results: bool,
                ici_mode: str):
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

    # Common RI index axis
    sip_df["rho_idx"]   = sip_df["rho"] / rho_ref
    sip_df["logRhoIdx"] = sip_df["logRho"] - logRho_ref

    # ICI index depends on mode (for GOAL plot only)
    if ici_mode == "sig_over_sigref":
        sip_df["imag_idx"]   = sip_df["sigma_imag"] / sig_ref
        sip_df["logImagIdx"] = sip_df["logImag"] - logSig_ref
    else:
        sip_df["imag_idx"]   = sig_ref / sip_df["sigma_imag"]
        sip_df["logImagIdx"] = logSig_ref - sip_df["logImag"]

    sip_df = sip_df[np.isfinite(sip_df["logRhoIdx"]) & np.isfinite(sip_df["logImagIdx"])]

    # RAW fits (logRho and logImag vs logS)
    slope_r, intercept_r, r2_r, p_r, ss_res_lin_r, aic_lin_r = linear_fit(sip_df["logS"], sip_df["logRho"])
    slope_i, intercept_i, r2_i, p_i, ss_res_lin_i, aic_lin_i = linear_fit(sip_df["logS"], sip_df["logImag"])
    lin_rho = {"slope": slope_r, "intercept": intercept_r, "r2": r2_r, "p": p_r, "aic": aic_lin_r}
    lin_im  = {"slope": slope_i, "intercept": intercept_i, "r2": r2_i, "p": p_i, "aic": aic_lin_i}

    # Through-origin fits on INDEX (goal plot)
    a_ri,  r2_ri,  _ = linear_fit_through_origin(sip_df["logS"], sip_df["logRhoIdx"])
    a_ici, r2_ici, _ = linear_fit_through_origin(sip_df["logS"], sip_df["logImagIdx"])

    # Report positive exponents consistently:
    n_idx = -a_ri
    if ici_mode == "sig_over_sigref":
        # logImagIdx = log(sig/sigref) => slope is +p
        p_idx = +a_ici
    else:
        # logImagIdx = log(sigref/sig) => slope is -p
        p_idx = -a_ici

    # ---- v27 TARGET (Van Nuys style): raw-like linear fits WITH intercept ----
    # Lock to VN definition always (independent of --ici-mode):
    sip_df["logRI_vn"]  = sip_df["logRhoIdx"]                 # log10(rho/rho_ref)
    sip_df["logICI_vn"] = logSig_ref - sip_df["logImag"]      # log10(sig_ref/sig_imag) => downtrend

    slope_RI_vn,  int_RI_vn,  r2_RI_vn,  _, _, _ = linear_fit(sip_df["logS"], sip_df["logRI_vn"])
    slope_ICI_vn, int_ICI_vn, r2_ICI_vn, _, _, _ = linear_fit(sip_df["logS"], sip_df["logICI_vn"])

    n_vn = -slope_RI_vn
    p_vn = -slope_ICI_vn

    sip_df["logRI_vn_fit"]  = slope_RI_vn  * sip_df["logS"] + int_RI_vn
    sip_df["logICI_vn_fit"] = slope_ICI_vn * sip_df["logS"] + int_ICI_vn

    # Logistic on RAW logRho vs logS
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

    # Residual columns for diagnostics
    logS = sip_df["logS"].to_numpy(float)
    logR = sip_df["logRho"].to_numpy(float)
    lin_pred = slope_r*logS + intercept_r
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

    # Write CSV
    csv_path = results_dir / f"joined_{tag}Hz.csv"
    sip_df[[
        "measurement","freq_hz","S","rho","sigma_imag",
        "logS","logRho","logImag",

        "rho_idx","imag_idx",
        "logRhoIdx","logImagIdx",

        # v27 TARGET VN columns
        "logRI_vn","logICI_vn","logRI_vn_fit","logICI_vn_fit",

        # diagnostics columns
        "resid_lin_rho","resid_log_rho",
        "logRho_linfit","logRho_logfit",
    ]].to_csv(csv_path, index=False)

    report = results_dir / f"fit_report_{tag}Hz.txt"
    report.write_text(
        f"Run root: {run_root}\n"
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

    # Plot scripts + gnuplot
    write_plot_target(results_dir, csv_path.name, f"target_exponents_{tag}Hz.png", target_freq,
                      n_vn, r2_RI_vn, p_vn, r2_ICI_vn)

    write_plot_index(results_dir, csv_path.name, f"index_exponents_{tag}Hz.png", target_freq,
                     n_idx, r2_ri, p_idx, r2_ici)

    write_plot_raw(results_dir, csv_path.name, f"exponents_raw_{tag}Hz.png", target_freq,
                   lin_rho, lin_im, log_fit_acc)

    write_plot_logistic_only(results_dir, csv_path.name, f"exponents_logistic_only_{tag}Hz.png", target_freq, log_fit_acc)

    write_plot_logistic_diag(results_dir, csv_path.name, f"exponents_logistic_diagnostics_{tag}Hz.png",
                             target_freq, has_log=(log_fit_acc is not None))

    if not dry_run:
        _run_gnuplot("plot_target_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_index_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_raw_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
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
        "ici_mode": ici_mode,
        "n_vn": float(n_vn) if np.isfinite(n_vn) else n_vn,
        "p_vn": float(p_vn) if np.isfinite(p_vn) else p_vn,
        "n_idx": float(n_idx) if np.isfinite(n_idx) else n_idx,
        "p_idx": float(p_idx) if np.isfinite(p_idx) else p_idx,
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
    ap.add_argument("--clean-results", action="store_true",
                    help="Delete files inside each run's Results/ before writing new outputs.")
    ap.add_argument("--ici-mode", type=str, default="sigref_over_sig",
                    choices=["sigref_over_sig","sig_over_sigref"],
                    help="GOAL/INDEX ICI definition: sigref_over_sig (default, ICI=sig_ref/sig) "
                         "or sig_over_sigref (ICI=sig/sig_ref). TARGET is always Van Nuys style.")
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
                run, args.freq, args.sample, args.dry_run,
                args.min_points, args.max_df,
                do_logistic=(not args.no_logistic),
                logistic_min_points=args.logistic_min_points,
                logistic_fast=args.logistic_fast,
                logistic_max_seconds=args.logistic_max_seconds,
                gnuplot_timeout=args.gnuplot_timeout,
                clean_results=args.clean_results,
                ici_mode=args.ici_mode
            )
            if res is None:
                continue
            rows.append(res)
            dt = time.time() - ti
            print(f"[{i}/{len(roots)}] OK points={res.get('points',0)} "
                  f"n_vn={res.get('n_vn',np.nan):.3g} p_vn={res.get('p_vn',np.nan):.3g} "
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
