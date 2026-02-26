#!/usr/bin/env python3
"""
auto_sip_pipeline_v18_noscipy.py

Fix vs v17:
- Compute INDEX columns using reference at highest saturation:
    logRhoIdx  = logRho  - logRho_ref_at_maxS
    logImagIdx = logImag - logImag_ref_at_maxS
  so index is 0 at S~1 (logS~0).
- Fit INDEX vs logS (default: through origin for physical consistency).
- Plot exponent figures in the "blue/red index" style (matches your goal).
- Keep logistic diagnostics infrastructure (still optional), but logistic is not used for index by default.
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

# Logistic constraints (avoid flat fits / overflow)
K_MIN = 0.05
K_MAX = 50.0

# Logistic acceptance settings
LOGISTIC_R2_EPS = 0.005
LOGISTIC_REQUIRE_BETTER_AIC = True

# ----------- small utils -----------
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

def normal_cdf(z):
    return 0.5*(1.0 + math.erf(z/math.sqrt(2.0)))

# ---------------- AIC ----------------
def aic_gaussian(ss_res, n, k_params):
    if n <= 0 or k_params <= 0 or not np.isfinite(ss_res) or ss_res <= 0:
        return np.nan
    return float(n * math.log(ss_res / n) + 2.0 * k_params)

# ---------------- Fits ----------------
def linear_fit(x, y):
    """OLS with intercept."""
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
    """OLS constrained through origin: y = a*x."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 2:
        return np.nan, np.nan, np.nan, np.nan
    denom = float(np.sum(x*x))
    if denom <= 0:
        return np.nan, np.nan, np.nan, np.nan
    slope = float(np.sum(x*y) / denom)
    yhat = slope*x
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    # F-test not strictly correct for constrained fit; use same helper as approximation with df_model=1
    pval = f_test_pvalue(ss_res, ss_tot, df_model=1, n=len(x))
    return slope, r2, pval, ss_res

# ---------------- Logistic (kept for diagnostics; acts on raw logRho vs logS) ----------------
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

    t0 = time.time()

    if fast:
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

        pvals = {"L": np.nan, "k": np.nan, "x0": np.nan, "c": np.nan}
        pvals_ok = False
        try:
            params = np.array([Lb, kb, x0b, cb], float)
            eps = 1e-5
            J = np.zeros((n, 4), float)
            for j in range(4):
                step = eps * max(1.0, abs(params[j]))
                dp = np.zeros(4); dp[j] = step
                y1 = logistic4(x, *(params + dp))
                y0 = logistic4(x, *(params - dp))
                J[:, j] = (y1 - y0) / (2.0 * step)

            dof = n - 4
            if dof > 0 and ss_res >= 0:
                sigma2 = ss_res / dof
                cov = sigma2 * np.linalg.pinv(J.T @ J)
                diag = np.diag(cov)
                if np.all(np.isfinite(diag)) and np.all(diag >= 0):
                    se = np.sqrt(diag)
                    if np.all(se > 0) and np.all(np.isfinite(se)):
                        t = params / se
                        for name, tj in zip(["L", "k", "x0", "c"], t):
                            pvals[name] = 2.0 * (1.0 - normal_cdf(abs(float(tj))))
                        pvals_ok = True
        except Exception:
            pvals_ok = False

        return {
            "L": float(Lb), "k": float(kb), "x0": float(x0b), "c": float(cb),
            "r2": float(r2),
            "p_model": float(p_model) if np.isfinite(p_model) else p_model,
            "aic": float(aic) if np.isfinite(aic) else aic,
            "pvals": pvals, "pvals_ok": pvals_ok,
            "ss_res": float(ss_res), "ss_tot": float(ss_tot)
        }, "ok"

    return None, "fail"

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

        out.append({"measurement": str(sh).strip(),
                    "freq_hz": f_found,
                    "rho": float(rho),
                    "sigma_imag": float(sim)})
    return pd.DataFrame(out)

# ---------------- Plot writing helpers ----------------
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

def write_gnuplot_index_style(results_dir: Path, csv_name: str, out_png: str,
                             target_freq: float,
                             slope_r_idx, r2_r_idx,
                             slope_i_idx, r2_i_idx):
    """
    Plots:
      y1 = log10(rho/rho_ref)  (>=0)
      y2 = log10(sig/sig_ref)  (<=0 typically)
    with lines forced through origin (index definition).
    """
    label_freq = f"{target_freq:g}"
    n_exp = -slope_r_idx
    p_exp = slope_i_idx

    gp = f"""reset
set datafile separator ','
set term pngcairo size 1200,1600 enhanced font 'Arial,28'
set output '{out_png}'
set multiplot layout 2,1

unset key
set grid
set xlabel 'log10(Saturation)'

# styles (match your goal)
set style line 1 lc rgb "#0000FF" lw 3 dt 2
set style line 2 lc rgb "#FF0000" lw 3 dt 2
set style line 3 lc rgb "#0000FF" pt 7 ps 1.8
set style line 4 lc rgb "#FF0000" pt 7 ps 1.8

# ---------------- Resistivity Index ----------------
set title 'Resistivity Index vs Saturation (log-log) at {label_freq} Hz'
set ylabel 'log10(Resistivity Index)'

a1 = {slope_r_idx:.12g}
f1(x) = a1*x

plot '{csv_name}' using 5:8 with points ls 3, \\
     f1(x) with lines ls 1

set label 1 sprintf('{label_freq}Hz: n=%.2f, R^2=%.3f', -a1, {r2_r_idx:.12g}) at graph 0.55,0.85 front

# ---------------- Imaginary Conductivity Index ----------------
set title 'Imaginary Conductivity Index vs Saturation (log-log) at {label_freq} Hz'
set ylabel 'log10(Imaginary Conductivity Index)'

a2 = {slope_i_idx:.12g}
f2(x) = a2*x

plot '{csv_name}' using 5:9 with points ls 4, \\
     f2(x) with lines ls 2

set label 2 sprintf('{label_freq}Hz: p=%.2f, R^2=%.3f', a2, {r2_i_idx:.12g}) at graph 0.55,0.85 front

unset multiplot
"""
    (results_dir / "plot_exponents.gp").write_text(gp)

def write_gnuplot_logistic_diagnostics(results_dir: Path, csv_name: str, out_png: str,
                                      target_freq: float, lin_rho, log_fit_accepted):
    """
    Diagnostics stays on RAW logRho vs logS (not index), because logistic is not meaningful on index (it is forced to 0 at S=1).
    """
    label_freq = f"{target_freq:g}"
    a_r = lin_rho["slope"]; b_r = lin_rho["intercept"]
    has_log = (log_fit_accepted is not None)

    log_defs = ""
    if has_log:
        L = log_fit_accepted["L"]; k = log_fit_accepted["k"]; x0 = log_fit_accepted["x0"]; c = log_fit_accepted["c"]
        log_defs = (
            f"L1 = {L:.12g}\n"
            f"k1 = {k:.12g}\n"
            f"x01= {x0:.12g}\n"
            f"c1 = {c:.12g}\n"
            f"fL(x) = c1 + L1/(1+exp(-k1*(x-x01)))\n"
        )

    panel1 = (
        f"set xrange [*:*]\nset yrange [*:*]\n"
        f"plot '{csv_name}' using 5:6 with points pt 7 ps 1.6 title 'data', \\\n"
        f"     f1(x) with lines dt 2 lw 3 title 'linear'"
        + ("" if not has_log else ", \\\n     fL(x) with lines dt 1 lw 4 title 'logistic (accepted)'")
        + "\n"
    )

    panel2 = (
        "set xrange [*:*]\nset yrange [*:*]\n"
        f"plot '{csv_name}' using 5:($6 - f1($5)) with points pt 7 ps 1.4 title 'residuals'\n"
    )

    if has_log:
        panel3 = (
            "set xrange [*:*]\nset yrange [*:*]\n"
            f"plot '{csv_name}' using 5:($6 - fL($5)) with points pt 7 ps 1.4 title 'residuals'\n"
        )
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
set term pngcairo size 1200,1800 enhanced font 'Arial,24'
set output '{out_png}'
set multiplot layout 3,1

a1 = {a_r:.12g}
b1 = {b_r:.12g}
f1(x) = a1*x + b1
{log_defs}

set title 'Diagnostics: logRho vs logS at {label_freq} Hz'
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

# ---------------- discovery / execution ----------------
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
                fit_through_origin: bool):
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")
    if sip_xlsx is None or log_xlsx is None:
        return None

    results_dir.mkdir(parents=True, exist_ok=True)

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

    sip_df = sip_df[np.isfinite(sip_df["logS"]) & np.isfinite(sip_df["logRho"]) & np.isfinite(sip_df["logImag"])].copy()

    if len(sip_df) < min_points:
        # still write whatever joined we have (useful debug)
        csv_path = results_dir / f"joined_{tag}Hz.csv"
        sip_df[["measurement","freq_hz","S","rho","logS","logRho","logImag"]].to_csv(csv_path, index=False)
        return {"run": str(run_root), "status": "FAIL",
                "reason": f"Not enough matched points (need >= {min_points})",
                "points": int(len(sip_df)), "logbook_sheet": used_sheet}

    # --- INDEX reference at max S (closest to saturation) ---
    ref_idx = int(np.nanargmax(sip_df["S"].to_numpy(float)))
    ref_logRho  = float(sip_df.iloc[ref_idx]["logRho"])
    ref_logImag = float(sip_df.iloc[ref_idx]["logImag"])

    sip_df["logRhoIdx"]  = sip_df["logRho"]  - ref_logRho
    sip_df["logImagIdx"] = sip_df["logImag"] - ref_logImag

    # write joined CSV including index columns
    csv_path = results_dir / f"joined_{tag}Hz.csv"
    sip_df[["measurement","freq_hz","S","rho","logS","logRho","logImag","logRhoIdx","logImagIdx"]].to_csv(csv_path, index=False)

    # ---- Fits on INDEX ----
    if fit_through_origin:
        slope_r, r2_r, p_r, ss_res_r = linear_fit_through_origin(sip_df["logS"], sip_df["logRhoIdx"])
        slope_i, r2_i, p_i, ss_res_i = linear_fit_through_origin(sip_df["logS"], sip_df["logImagIdx"])
        intercept_r = 0.0
        intercept_i = 0.0
        aic_r = aic_gaussian(ss_res_r, len(sip_df), k_params=1)
        aic_i = aic_gaussian(ss_res_i, len(sip_df), k_params=1)
    else:
        slope_r, intercept_r, r2_r, p_r, ss_res_r, aic_r = linear_fit(sip_df["logS"], sip_df["logRhoIdx"])
        slope_i, intercept_i, r2_i, p_i, ss_res_i, aic_i = linear_fit(sip_df["logS"], sip_df["logImagIdx"])

    n_exp = -slope_r
    p_exp = slope_i

    # ---- Optional logistic on RAW (diagnostics only) ----
    log_fit = None
    log_status = "skipped"
    if do_logistic:
        if len(sip_df) >= logistic_min_points:
            log_fit, log_status = logistic_fit_noscipy(
                sip_df["logS"], sip_df["logRho"],
                fast=logistic_fast,
                max_seconds=(logistic_max_seconds if logistic_max_seconds and logistic_max_seconds > 0 else 0.0)
            )
        else:
            log_status = "not_enough_points"

    log_fit_acc = None
    log_accept_reason = "n/a"
    if log_fit is not None:
        # linear baseline for raw logRho (only for acceptance logic)
        slope_raw, intercept_raw, r2_raw, p_raw, ss_raw, aic_raw = linear_fit(sip_df["logS"], sip_df["logRho"])
        ok, reason = accept_logistic(r2_raw, aic_raw, log_fit)
        log_accept_reason = reason
        if ok:
            log_fit_acc = log_fit

    # ---- report ----
    report = results_dir / f"fit_report_{tag}Hz.txt"
    txt = (
        f"Run root: {run_root}\n"
        f"SIP file: {sip_xlsx}\n"
        f"Logbook : {log_xlsx}\n"
        f"Logbook sheet used: {used_sheet}\n"
        f"Target freq requested: {target_freq} Hz\n"
        f"Matched points used: {len(sip_df)}\n"
        f"Index reference point (max S): measurement={sip_df.iloc[ref_idx]['measurement']}, "
        f"S={sip_df.iloc[ref_idx]['S']:.6g}\n\n"
        f"--- INDEX fit (logRhoIdx vs logS) ---\n"
        f"logRhoIdx = a*logS + b\n"
        f"a = {slope_r:.6g}\n"
        f"b = {intercept_r:.6g}\n"
        f"n = {-slope_r:.6g}\n"
        f"R2 = {r2_r:.6g}\n"
        f"p(model) = {p_r:.6g}\n"
        f"SSE = {ss_res_r:.6g}\n"
        f"AIC = {aic_r:.6g}\n\n"
        f"--- INDEX fit (logImagIdx vs logS) ---\n"
        f"logImagIdx = p*logS + q\n"
        f"p = {slope_i:.6g}\n"
        f"q = {intercept_i:.6g}\n"
        f"R2 = {r2_i:.6g}\n"
        f"p(model) = {p_i:.6g}\n"
        f"SSE = {ss_res_i:.6g}\n"
        f"AIC = {aic_i:.6g}\n\n"
        f"(Derived exponents) n={n_exp:.6g}, p={p_exp:.6g}\n\n"
    )

    if log_fit is None:
        txt += f"Logistic (raw logRho): {log_status}\n"
    else:
        txt += (
            f"--- Logistic 4p fit on RAW logRho vs logS (diagnostics only) ---\n"
            f"status = {log_status}\n"
            f"L  = {log_fit['L']:.6g}\n"
            f"k  = {log_fit['k']:.6g}\n"
            f"x0 = {log_fit['x0']:.6g}\n"
            f"c  = {log_fit['c']:.6g}\n"
            f"R2 = {log_fit['r2']:.6g}\n"
            f"p(model) = {log_fit['p_model']:.6g}\n"
            f"SSE = {log_fit.get('ss_res', np.nan):.6g}\n"
            f"AIC = {log_fit.get('aic', np.nan):.6g}\n"
            f"ACCEPTED = {bool(log_fit_acc is not None)}\n"
            f"accept_reason = {log_accept_reason}\n"
        )

    report.write_text(txt)

    # ---- plots ----
    out_png = f"exponents_{tag}Hz.png"
    write_gnuplot_index_style(
        results_dir, csv_path.name, out_png, target_freq,
        slope_r, r2_r,
        slope_i, r2_i
    )

    out_png_diag = f"exponents_logistic_diagnostics_{tag}Hz.png"
    # linear for diagnostics uses raw logRho
    slope_raw, intercept_raw, r2_raw, p_raw, ss_raw, aic_raw = linear_fit(sip_df["logS"], sip_df["logRho"])
    lin_raw = {"slope": slope_raw, "intercept": intercept_raw, "r2": r2_raw, "p": p_raw, "aic": aic_raw}
    write_gnuplot_logistic_diagnostics(results_dir, csv_path.name, out_png_diag, target_freq, lin_raw, log_fit_acc)

    if not dry_run:
        _run_gnuplot("plot_exponents.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=True)
        _run_gnuplot("plot_logistic_diagnostics.gp", cwd=results_dir, timeout_s=gnuplot_timeout, fatal=False)

    return {"run": str(run_root), "status": "OK", "reason": "",
            "points": int(len(sip_df)), "logbook_sheet": used_sheet,
            "n": float(n_exp) if np.isfinite(n_exp) else n_exp,
            "p": float(p_exp) if np.isfinite(p_exp) else p_exp,
            "r2_rho_idx": float(r2_r) if np.isfinite(r2_r) else r2_r,
            "r2_imag_idx": float(r2_i) if np.isfinite(r2_i) else r2_i,
            "logistic_attempted": bool(log_fit is not None),
            "logistic_accepted": bool(log_fit_acc is not None),
            "logistic_status": log_status,
            "logistic_reason": log_accept_reason}

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

    # NEW: for index fits
    ap.add_argument("--fit-intercept", action="store_true",
                    help="Fit index with intercept (default fits through origin).")

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
                fit_through_origin=(not args.fit_intercept)
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
