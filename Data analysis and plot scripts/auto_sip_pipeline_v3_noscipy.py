#!/usr/bin/env python3
import os, sys, math, argparse, subprocess
from pathlib import Path

import numpy as np
import pandas as pd

# -----------------------------
# Config defaults
# -----------------------------
ANALYSIS_DIR = "Analysis"
RAW_DIR      = "SIP_Raw_Data"
RESULTS_DIR  = "Results"

ANALYSIS_XLSX_HINT = "output"
LOGBOOK_XLSX_HINT  = "Log Book"

TARGET_FREQ_DEFAULT = 0.01

LOGBOOK_ID_COL = "Measurements"
LOGBOOK_SAT_COL = "Degree Saturation M2"
LOGBOOK_SAMPLE_COL = "Sample"

# SIP output typical headers (we also fallback by position)
SIP_FREQ_CANDS = ["Frequency", "Frequency (Hz)", "Freq", "Hz"]
SIP_RHO_CANDS  = ["Resistivity", "Resistivity (Ω·m)", "Resistivity (Ohm m)", "rho"]
SIP_IMAG_CANDS = ["Imaginary Conductivity", "Imaginary Conductivity (µS/cm)", "Imag", "sigma_imag"]

def find_first_matching_file(folder: Path, hint: str, ext=".xlsx"):
    if not folder.exists(): return None
    for p in folder.rglob(f"*{ext}"):
        if hint.lower() in p.name.lower():
            return p
    return None

def pick_col(df: pd.DataFrame, cands):
    low = {c.lower(): c for c in df.columns}
    for k in cands:
        if k.lower() in low:
            return low[k.lower()]
    return ""

def to_float(x):
    try:
        if isinstance(x, str):
            xs = x.strip().replace("%","")
            return float(xs)
        return float(x)
    except Exception:
        return float("nan")

def safe_log10(x):
    x = float(x)
    if not np.isfinite(x) or x <= 0: return np.nan
    return math.log10(x)

def linear_fit(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 2:
        return np.nan, np.nan, np.nan, np.nan
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    yhat = slope*x + intercept
    ss_res = np.sum((y - yhat)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    # F-test p-value vs intercept-only
    # df1 = 1, df2 = n-2
    pval = f_test_pvalue(ss_res, ss_tot, df_model=1, n=len(x))
    return slope, intercept, r2, pval

def logistic4(x, L, k, x0, c):
    # y = c + L/(1+exp(-k(x-x0)))
    return c + L / (1.0 + np.exp(-k*(x - x0)))

def sse(y, yhat):
    r = y - yhat
    return float(np.sum(r*r))

def f_test_pvalue(ss_res, ss_tot, df_model, n):
    # Compare to intercept-only model: ss_tot = SSE0
    # F = ((SSE0 - SSE1)/df_model) / (SSE1/(n-df_model-1))
    # df2 = n - df_model - 1
    df2 = n - df_model - 1
    if df2 <= 0 or not np.isfinite(ss_res) or not np.isfinite(ss_tot): 
        return np.nan
    num = (ss_tot - ss_res)/df_model
    den = ss_res/df2
    if den <= 0: 
        return np.nan
    F = num/den
    # p = 1 - CDF_F(F; df_model, df2)
    return f_dist_sf(F, df_model, df2)

# ---- F distribution survival function using incomplete beta (no SciPy) ----
def betacf(a,b,x, maxit=200, eps=3e-14):
    # Continued fraction for incomplete beta (Numerical Recipes)
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

def betai(a,b,x):
    # Regularized incomplete beta I_x(a,b)
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    ln_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a+b)
    bt = math.exp(math.log(x)*a + math.log(1.0-x)*b - ln_beta)
    if x < (a+1.0)/(a+b+2.0):
        return bt*betacf(a,b,x)/a
    else:
        return 1.0 - bt*betacf(b,a,1.0-x)/b

def f_dist_sf(F, d1, d2):
    # Survival function for F(d1,d2)
    # CDF = I_{d1*F/(d1*F+d2)}(d1/2, d2/2)
    if F < 0: return 1.0
    x = (d1*F)/(d1*F + d2)
    cdf = betai(d1/2.0, d2/2.0, x)
    return max(0.0, 1.0 - cdf)

def logistic_fit_noscipy(x, y):
    """
    Grid search + local refinement for (L,k,x0,c).
    Returns params, R2, p_model, per-param pvals (approx) and SSE.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    n = len(x)
    if n < 5:
        return None

    # initial ranges
    y_min, y_max = float(np.min(y)), float(np.max(y))
    x_min, x_max = float(np.min(x)), float(np.max(x))

    # heuristic initial guesses
    L0 = (y_max - y_min) if (y_max>y_min) else 1.0
    c0 = y_min
    x00 = float(np.median(x))
    k0 = 2.0

    # coarse grid
    L_grid  = np.linspace(0.3*L0, 2.0*L0, 18)
    c_grid  = np.linspace(y_min - 0.5*L0, y_min + 0.5*L0, 14)
    x0_grid = np.linspace(x_min, x_max, 16)
    k_grid  = np.concatenate([np.linspace(0.2, 2.0, 10), np.linspace(2.5, 10.0, 8)])

    best = None
    for L in L_grid:
        for c in c_grid:
            for x0 in x0_grid:
                for k in k_grid:
                    yhat = logistic4(x, L, k, x0, c)
                    s = sse(y, yhat)
                    if (best is None) or (s < best[0]):
                        best = (s, L, k, x0, c)

    _, Lb, kb, x0b, cb = best

    # local refinement: coordinate descent
    stepL = 0.15*max(abs(Lb), 1e-6)
    stepk = 0.15*max(abs(kb), 1e-6)
    stepx0= 0.10*max(abs(x_max-x_min), 1e-6)
    stepc = 0.15*max(abs(Lb), 1e-6)

    def try_update(L,k,x0,c, stepL,stepk,stepx0,stepc):
        base = sse(y, logistic4(x, L,k,x0,c))
        improved = True
        while improved:
            improved = False
            for dL in [0, -stepL, stepL]:
                for dk in [0, -stepk, stepk]:
                    for dx0 in [0, -stepx0, stepx0]:
                        for dc in [0, -stepc, stepc]:
                            if dL==dk==dx0==dc==0: continue
                            L2 = max(1e-9, L + dL)
                            k2 = max(1e-6, k + dk)
                            x02 = x0 + dx0
                            c2 = c + dc
                            s2 = sse(y, logistic4(x, L2,k2,x02,c2))
                            if s2 + 1e-12 < base:
                                base = s2
                                L,k,x0,c = L2,k2,x02,c2
                                improved = True
        return L,k,x0,c,base

    for _ in range(6):
        Lb,kb,x0b,cb, s_best = try_update(Lb,kb,x0b,cb, stepL,stepk,stepx0,stepc)
        stepL *= 0.5; stepk *= 0.5; stepx0 *= 0.5; stepc *= 0.5

    yhat = logistic4(x, Lb, kb, x0b, cb)
    ss_res = sse(y, yhat)
    ss_tot = sse(y, np.full_like(y, np.mean(y)))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else np.nan
    p_model = f_test_pvalue(ss_res, ss_tot, df_model=3, n=n)  # logistic vs intercept-only

    # Approx parameter p-values from Jacobian (finite differences)
    pvals = {"L": np.nan, "k": np.nan, "x0": np.nan, "c": np.nan}
    try:
        params = np.array([Lb,kb,x0b,cb], float)
        eps = 1e-5
        J = np.zeros((n,4), float)
        for j in range(4):
            dp = np.zeros(4); dp[j] = eps*max(1.0, abs(params[j]))
            y1 = logistic4(x, *(params + dp))
            y0 = logistic4(x, *(params - dp))
            J[:,j] = (y1 - y0) / (2.0*dp[j])
        # covariance ~ sigma^2 * (J^T J)^-1
        dof = n - 4
        if dof > 0:
            sigma2 = ss_res / dof
            JTJ = J.T @ J
            cov = sigma2 * np.linalg.inv(JTJ)
            se = np.sqrt(np.diag(cov))
            t = params / se
            # two-sided p from t approx normal for small n -> use normal approx
            # (with very small n, this is approximate; still informative)
            for name, tj in zip(["L","k","x0","c"], t):
                pvals[name] = 2.0*(1.0 - normal_cdf(abs(tj)))
    except Exception:
        pass

    return {
        "L": Lb, "k": kb, "x0": x0b, "c": cb,
        "r2": r2, "p_model": p_model, "pvals": pvals, "sse": ss_res
    }

def normal_cdf(z):
    # CDF of standard normal
    return 0.5*(1.0 + math.erf(z/math.sqrt(2.0)))

def read_logbook(log_xlsx: Path, sample_filter: str|None):
    df = pd.read_excel(log_xlsx)
    if sample_filter and (LOGBOOK_SAMPLE_COL in df.columns):
        df = df[df[LOGBOOK_SAMPLE_COL].astype(str).str.contains(sample_filter, case=False, na=False)]
    if (LOGBOOK_ID_COL not in df.columns) or (LOGBOOK_SAT_COL not in df.columns):
        raise RuntimeError(f"Logbook missing columns. Need '{LOGBOOK_ID_COL}' and '{LOGBOOK_SAT_COL}'. Got: {list(df.columns)}")
    sat_map = {}
    for _, r in df.iterrows():
        mid = str(r[LOGBOOK_ID_COL]).strip()
        S = to_float(r[LOGBOOK_SAT_COL])
        if np.isfinite(S):
            if S > 1.0 and S <= 100.0:  # percent
                S = S/100.0
        if mid and np.isfinite(S):
            sat_map[mid] = float(S)
    return sat_map

def read_sip_output(sip_xlsx: Path, target_freq: float):
    xls = pd.ExcelFile(sip_xlsx)
    out = []
    for sh in xls.sheet_names:
        df = pd.read_excel(sip_xlsx, sheet_name=sh)
        freq_col = pick_col(df, SIP_FREQ_CANDS) or df.columns[0]
        rho_col  = pick_col(df, SIP_RHO_CANDS)
        imag_col = pick_col(df, SIP_IMAG_CANDS)

        if (not rho_col) or (not imag_col):
            # fallback to your known format: freq, imag, phase, real, resistivity
            if len(df.columns) >= 5:
                imag_col = df.columns[1]
                rho_col  = df.columns[4]
            else:
                continue

        freqs = df[freq_col].apply(to_float).to_numpy(float)
        idx = int(np.nanargmin(np.abs(freqs - target_freq)))
        f_found = freqs[idx]
        rho = to_float(df.loc[df.index[idx], rho_col])
        sim = to_float(df.loc[df.index[idx], imag_col])
        out.append({
            "measurement": str(sh).strip(),
            "freq_hz": float(f_found),
            "rho": float(rho),
            "sigma_imag": float(sim)
        })
    return pd.DataFrame(out)

def write_gnuplot(results_dir: Path, csv_name: str, out_png: str,
                  n_lin, r2_lin, p_lin,
                  p_log, r2_log, pval_log, L,k,x0,c,
                  label_freq: str):

    gp = f"""\
reset
set datafile separator ','
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'

set multiplot layout 2,1

# Panel 1: Resistivity
set title 'Van Nuys S1: Resistivity vs Saturation at {label_freq} Hz'
set grid
set key top right
set xlabel 'log(Saturation)'
set ylabel 'log(Resistivity)'

# Linear fit line: y = a*x + b, with n = -a
a1 = {-n_lin:.12g}
b1 = 0
f1(x) = a1*x + b1

# Logistic curve parameters on log-log space
L1 = {L:.12g}
k1 = {k:.12g}
x01= {x0:.12g}
c1 = {c:.12g}
fL(x) = c1 + L1/(1+exp(-k1*(x-x01)))

plot '{csv_name}' using 5:6 with points pt 7 ps 1.7 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('linear: n=%.2f, R^2=%.3f, p=%.3g', -a1, {r2_lin:.12g}, {p_lin:.12g}), \\
     fL(x) with lines dt 1 lw 4 title sprintf('logistic: L=%.2f k=%.2f x0=%.2f c=%.2f, R^2=%.3f, p=%.3g', L1,k1,x01,c1, {r2_log:.12g}, {pval_log:.12g})

# Panel 2: Imag cond
set title 'Van Nuys S1: Imag Cond vs Saturation at {label_freq} Hz'
set grid
set key top right
set xlabel 'log(Saturation)'
set ylabel 'log(Imaginary Conductivity)'

a2 = {p_log:.12g}
b2 = 0
f2(x) = a2*x + b2

plot '{csv_name}' using 5:7 with points pt 7 ps 1.7 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('linear: p=%.2f, R^2=%.3f', a2, {r2_log:.12g})

unset multiplot
"""
    (results_dir / "plot_exponents.gp").write_text(gp)

def process_run(run_root: Path, target_freq: float, sample_filter: str|None, dry_run: bool):
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir      = run_root / RAW_DIR
    results_dir  = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")

    if sip_xlsx is None or log_xlsx is None:
        return False

    results_dir.mkdir(parents=True, exist_ok=True)

    sat_map = read_logbook(log_xlsx, sample_filter)
    sip_df = read_sip_output(sip_xlsx, target_freq)

    if sip_df.empty:
        raise RuntimeError(f"No SIP sheets read from {sip_xlsx}")

    sip_df["S"] = sip_df["measurement"].map(sat_map).astype(float)

    sip_df["logS"] = sip_df["S"].apply(safe_log10)
    sip_df["logRho"] = sip_df["rho"].apply(safe_log10)
    sip_df["logImag"] = sip_df["sigma_imag"].apply(safe_log10)

    # Drop rows without saturation
    sip_df = sip_df[np.isfinite(sip_df["logS"]) & np.isfinite(sip_df["logRho"]) & np.isfinite(sip_df["logImag"])]

    if len(sip_df) < 4:
        raise RuntimeError(f"Not enough matched points after join (need >=4). Got {len(sip_df)}. Check logbook mapping Measurements <-> sheet names.")

    csv_path = results_dir / f"S1_joined_{str(target_freq).replace('.','p')}Hz.csv"
    sip_df[["measurement","freq_hz","S","rho","logS","logRho","logImag"]].to_csv(csv_path, index=False)

    # Linear fits
    slope_r, intercept_r, r2_r, p_r = linear_fit(sip_df["logS"], sip_df["logRho"])
    n = -slope_r
    slope_i, intercept_i, r2_i, p_i = linear_fit(sip_df["logS"], sip_df["logImag"])

    # Logistic fit on resistivity (log-log space)
    fit = logistic_fit_noscipy(sip_df["logS"], sip_df["logRho"])
    if fit is None:
        raise RuntimeError("Not enough points for logistic fit (need >=5).")

    # Write report
    rep = results_dir / f"S1_fit_report_{str(target_freq).replace('.','p')}Hz.txt"
    rep.write_text(
        f"Run root: {run_root}\n"
        f"SIP file: {sip_xlsx}\n"
        f"Logbook : {log_xlsx}\n"
        f"Target freq requested: {target_freq} Hz\n"
        f"Points used: {len(sip_df)}\n\n"
        f"--- Linear fit (logRho vs logS): logRho = a*logS + b ---\n"
        f"a = {slope_r:.6g}\n"
        f"b = {intercept_r:.6g}\n"
        f"n = {-slope_r:.6g}\n"
        f"R2 = {r2_r:.6g}\n"
        f"p(model) = {p_r:.6g}\n\n"
        f"--- Linear fit (logImag vs logS): logImag = p*logS + b ---\n"
        f"p = {slope_i:.6g}\n"
        f"b = {intercept_i:.6g}\n"
        f"R2 = {r2_i:.6g}\n"
        f"p(model) = {p_i:.6g}\n\n"
        f"--- Logistic 4p fit on logRho vs logS: y = c + L/(1+exp(-k(x-x0))) ---\n"
        f"L  = {fit['L']:.6g}\n"
        f"k  = {fit['k']:.6g}\n"
        f"x0 = {fit['x0']:.6g}\n"
        f"c  = {fit['c']:.6g}\n"
        f"R2 = {fit['r2']:.6g}\n"
        f"p(model) = {fit['p_model']:.6g}\n"
        f"p(L)  ~ {fit['pvals']['L']:.6g}\n"
        f"p(k)  ~ {fit['pvals']['k']:.6g}\n"
        f"p(x0) ~ {fit['pvals']['x0']:.6g}\n"
        f"p(c)  ~ {fit['pvals']['c']:.6g}\n"
    )

    out_png = f"S1_exponents_{str(target_freq).replace('.','p')}Hz.png"
    label_freq = f"{target_freq:g}"

    write_gnuplot(
        results_dir=results_dir,
        csv_name=csv_path.name,
        out_png=out_png,
        n_lin=n, r2_lin=r2_r, p_lin=p_r,
        p_log=slope_i, r2_log=r2_i, pval_log=fit["p_model"],
        L=fit["L"], k=fit["k"], x0=fit["x0"], c=fit["c"],
        label_freq=label_freq
    )

    if not dry_run:
        subprocess.run(["gnuplot", "plot_exponents.gp"], cwd=str(results_dir), check=True)

    print(f"[OK] {run_root} -> {results_dir/out_png}")
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("top", help="Top parent folder to scan recursively")
    ap.add_argument("--freq", type=float, default=TARGET_FREQ_DEFAULT, help="Target frequency in Hz (default 0.01)")
    ap.add_argument("--sample", type=str, default=None, help="Optional sample filter (substring match in logbook Sample column)")
    ap.add_argument("--dry-run", action="store_true", help="Do not run gnuplot (still writes Results files)")
    args = ap.parse_args()

    base = Path(args.top).expanduser().resolve()
    if not base.exists():
        print(f"Not found: {base}")
        sys.exit(2)

    count = 0
    for p in base.rglob("*"):
        if p.is_dir() and (p/ANALYSIS_DIR).is_dir() and (p/RAW_DIR).is_dir():
            try:
                if process_run(p, args.freq, args.sample, args.dry_run):
                    count += 1
            except Exception as e:
                print(f"[FAIL] {p}: {e}")

    print(f"Done. Processed {count} run folder(s).")

if __name__ == "__main__":
    main()
