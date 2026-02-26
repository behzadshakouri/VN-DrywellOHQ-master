#!/usr/bin/env python3
"""
Fully automatic SIP + saturation pipeline (recursive).

For each "run root" directory that contains:
  - Analysis/      (contains an SIP output .xlsx with sheets M1, M2, ...)
  - SIP_Raw_Data/  (contains a logbook .xlsx with saturation information)

This script will:
  1) Create a sibling folder: Results/   (does NOT touch existing folders/files)
  2) Read degree of saturation S for each measurement (M1, M2, ...) from the logbook
  3) Read resistivity rho and imaginary conductivity sigma_imag at a target frequency from the SIP output workbook
  4) Join them on measurement ID
  5) Fit BOTH:
       - log-log linear model   (Archie-style): log(y) = a*log(S) + b
       - logistic curve on logs: y = c + L / (1 + exp(-k*(x - x0)))
         where x = log(S) and y = log(rho) or log(sigma_imag)
     and report R², model p-value, and parameters (L, k, x0, c).
  6) Write:
       - joined CSV
       - a small TXT report (fits + params)
       - a gnuplot script
       - the output PNG
     into Results/

Usage:
  python3 auto_sip_pipeline.py /path/to/parent --freq 0.01
  python3 auto_sip_pipeline.py /path/to/parent --freq 0.01 --sample "Van Nuys S-1"
  python3 auto_sip_pipeline.py /path/to/parent --freq 0.01 --dry-run
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats


# -----------------------------
# Defaults tailored to your uploaded files
# -----------------------------
ANALYSIS_DIR = "Analysis"
RAW_DIR = "SIP_Raw_Data"
RESULTS_DIR = "Results"

# file hints
ANALYSIS_XLSX_HINT = "output"
LOGBOOK_XLSX_HINT = "Log Book"

# logbook sheet + columns (from your header row)
LOGBOOK_SHEET = "Saturation"
LOGBOOK_ID_COL = "Measurements"
LOGBOOK_SAMPLE_COL = "Sample"
LOGBOOK_SAT_COL = "Degree Saturation M2"

# SIP output columns (from your uploaded workbook)
SIP_FREQ_COL = "Frequency (Hz)"
SIP_RHO_COL = "Resistivity (Ω·m)"
SIP_IMAG_COL = "Imaginary_Conductivity (µS/cm)"


def find_first_matching_file(folder: Path, must_contain: str, ext: str = ".xlsx") -> Path | None:
    if not folder.exists():
        return None
    for p in folder.rglob(f"*{ext}"):
        if must_contain.lower() in p.name.lower():
            return p
    return None


def safe_log10(v: float) -> float:
    try:
        v = float(v)
    except Exception:
        return float("nan")
    if not math.isfinite(v) or v <= 0:
        return float("nan")
    return math.log10(v)


def logistic4(x, L, k, x0, c):
    # y = c + L / (1 + exp(-k*(x-x0)))
    return c + L / (1.0 + np.exp(-k * (x - x0)))


def fit_logistic(x, y):
    """
    Fit logistic4 to (x,y). Returns dict with params, R2, p_model, param_pvals.
    Uses an approximate F-test vs intercept-only, with df based on parameter counts.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    out = {
        "L": float("nan"), "k": float("nan"), "x0": float("nan"), "c": float("nan"),
        "R2": float("nan"), "p_model": float("nan"),
        "p_L": float("nan"), "p_k": float("nan"), "p_x0": float("nan"), "p_c": float("nan"),
        "n": int(len(x))
    }
    if len(x) < 4:
        return out

    # initial guesses
    y_min = float(np.nanmin(y))
    y_max = float(np.nanmax(y))
    c0 = y_min
    L0 = max(1e-6, y_max - y_min)
    x0_0 = float(np.nanmedian(x))
    k0 = 1.0

    # bounds: keep L positive/negative allowed? Usually positive; allow either by wide bounds.
    # Use moderate bounds to avoid overflow.
    bounds = ([-1e6, -1e3, -10.0, -1e6],
              [ 1e6,  1e3,  10.0,  1e6])

    try:
        popt, pcov = curve_fit(
            logistic4, x, y,
            p0=[L0, k0, x0_0, c0],
            bounds=bounds,
            maxfev=20000
        )
        L, k, x0, c = [float(v) for v in popt]
        out.update({"L": L, "k": k, "x0": x0, "c": c})
    except Exception:
        return out

    # predictions + R2
    yhat = logistic4(x, out["L"], out["k"], out["x0"], out["c"])
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    out["R2"] = R2

    # model p-value via approximate F-test vs intercept-only
    # full model: 4 params, null model: 1 param (mean)
    n = len(x)
    p_full = 4
    p_null = 1
    df1 = p_full - p_null
    df2 = n - p_full
    if df2 > 0 and math.isfinite(ss_tot) and ss_tot > 0 and math.isfinite(ss_res):
        ss_null = ss_tot  # residuals of intercept-only model
        F = ((ss_null - ss_res) / df1) / (ss_res / df2) if ss_res > 0 else float("inf")
        p_model = 1.0 - stats.f.cdf(F, df1, df2) if math.isfinite(F) else float("nan")
        out["p_model"] = float(p_model)

    # parameter p-values (approx t-test)
    try:
        se = np.sqrt(np.diag(pcov))
        df = max(1, n - p_full)
        tvals = popt / se
        pvals = 2.0 * (1.0 - stats.t.cdf(np.abs(tvals), df))
        out["p_L"], out["p_k"], out["p_x0"], out["p_c"] = [float(v) for v in pvals]
    except Exception:
        pass

    return out


def fit_linear(x, y):
    """
    y = a*x + b. Returns a, b, R2, p_model.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    out = {"a": float("nan"), "b": float("nan"), "R2": float("nan"), "p_model": float("nan"), "n": int(len(x))}
    if len(x) < 2:
        return out

    a, b = np.polyfit(x, y, 1)
    yhat = a * x + b
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    out.update({"a": float(a), "b": float(b), "R2": float(R2)})

    # p-value for slope via Pearson correlation test (equivalent for simple linear regression)
    if len(x) >= 3:
        r, p = stats.pearsonr(x, y)
        out["p_model"] = float(p)
    return out


def freq_to_tag(freq: float) -> str:
    # 0.01 -> 0p01Hz, 10 -> 10Hz
    s = f"{freq:g}"
    s = s.replace(".", "p")
    return f"{s}Hz"


def write_gnuplot(out_dir: Path, csv_name: str, out_png: str, freq: float,
                 lin_r, logi_r, lin_i, logi_i):
    """
    Two panels:
      - log(rho) vs log(S) with linear + logistic fits
      - log(imag) vs log(S) with linear + logistic fits
    """
    # Embed parameters for plotting
    gp = f"""\
reset
set datafile separator ','
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'

set multiplot layout 2,1

# Common
set grid
set key top right
set xlabel 'log(Saturation)'

# --------------------
# Panel 1: Resistivity
# --------------------
set title 'Van Nuys S1: log(Resistivity) vs log(S) at {freq:g} Hz'
set ylabel 'log(Resistivity)'

# Linear fit: y = a*x + b  (Archie slope => n = -a)
a1 = {lin_r["a"]:.12g}
b1 = {lin_r["b"]:.12g}
f1(x) = a1*x + b1

# Logistic fit: y = c + L/(1+exp(-k*(x-x0)))
L1 = {logi_r["L"]:.12g}
k1 = {logi_r["k"]:.12g}
x01 = {logi_r["x0"]:.12g}
c1 = {logi_r["c"]:.12g}
g1(x) = c1 + L1/(1+exp(-k1*(x-x01)))

plot '{csv_name}' using 5:6 with points pt 7 ps 1.6 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('linear: n=%.2f, R^2=%.3f, p=%.3g', -a1, {lin_r["R2"]:.6g}, {lin_r["p_model"]:.6g}), \\
     g1(x) with lines dt 1 lw 3 title sprintf('logistic: R^2=%.3f, p=%.3g (L=%.2f,k=%.2f,x0=%.2f,c=%.2f)', {logi_r["R2"]:.6g}, {logi_r["p_model"]:.6g}, L1, k1, x01, c1)

# --------------------
# Panel 2: Imag cond
# --------------------
set title 'Van Nuys S1: log(Imag. Conductivity) vs log(S) at {freq:g} Hz'
set ylabel 'log(Imag. Conductivity)'

a2 = {lin_i["a"]:.12g}
b2 = {lin_i["b"]:.12g}
f2(x) = a2*x + b2

L2 = {logi_i["L"]:.12g}
k2 = {logi_i["k"]:.12g}
x02 = {logi_i["x0"]:.12g}
c2 = {logi_i["c"]:.12g}
g2(x) = c2 + L2/(1+exp(-k2*(x-x02)))

plot '{csv_name}' using 5:7 with points pt 7 ps 1.6 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('linear: p=%.2f, R^2=%.3f, p=%.3g', a2, {lin_i["R2"]:.6g}, {lin_i["p_model"]:.6g}), \\
     g2(x) with lines dt 1 lw 3 title sprintf('logistic: R^2=%.3f, p=%.3g (L=%.2f,k=%.2f,x0=%.2f,c=%.2f)', {logi_i["R2"]:.6g}, {logi_i["p_model"]:.6g}, L2, k2, x02, c2)

unset multiplot
"""
    (out_dir / "plot_exponents.gp").write_text(gp)


def read_logbook_saturation(log_path: Path, sample_filter: str | None):
    """
    Returns dict: measurement_id -> S (0..1)
    Handles:
      - forward-fill Sample column when blank
      - Degree Saturation in percent (0..100) or fraction (0..1)
    """
    df = pd.read_excel(log_path, sheet_name=LOGBOOK_SHEET)

    # forward fill sample names if needed
    if LOGBOOK_SAMPLE_COL in df.columns:
        df[LOGBOOK_SAMPLE_COL] = df[LOGBOOK_SAMPLE_COL].ffill()

    if sample_filter:
        # case-insensitive substring match
        mask = df[LOGBOOK_SAMPLE_COL].astype(str).str.contains(sample_filter, case=False, na=False)
        df = df[mask].copy()

    if LOGBOOK_ID_COL not in df.columns or LOGBOOK_SAT_COL not in df.columns:
        raise RuntimeError(
            f"Logbook missing expected columns. Need '{LOGBOOK_ID_COL}' and '{LOGBOOK_SAT_COL}'. "
            f"Have: {list(df.columns)}"
        )

    sat_map = {}
    for _, r in df.iterrows():
        mid = str(r[LOGBOOK_ID_COL]).strip()
        if not mid or mid.lower() == "nan":
            continue
        s = r[LOGBOOK_SAT_COL]
        try:
            s = float(s)
        except Exception:
            continue
        if not math.isfinite(s):
            continue
        # convert percent to fraction if needed
        if s > 1.5:
            s = s / 100.0
        sat_map[mid] = s

    return sat_map


def extract_sip_at_freq(sip_path: Path, target_freq: float):
    """
    Returns list of dict rows for each sheet measurement at target freq.
    """
    xls = pd.ExcelFile(sip_path)
    rows = []

    for sh in xls.sheet_names:
        df = pd.read_excel(sip_path, sheet_name=sh)

        if SIP_FREQ_COL not in df.columns:
            raise RuntimeError(f"SIP sheet '{sh}' missing '{SIP_FREQ_COL}'. Have {list(df.columns)}")
        if SIP_RHO_COL not in df.columns:
            raise RuntimeError(f"SIP sheet '{sh}' missing '{SIP_RHO_COL}'. Have {list(df.columns)}")
        if SIP_IMAG_COL not in df.columns:
            raise RuntimeError(f"SIP sheet '{sh}' missing '{SIP_IMAG_COL}'. Have {list(df.columns)}")

        freqs = pd.to_numeric(df[SIP_FREQ_COL], errors="coerce")
        # nearest row
        idx = (freqs - target_freq).abs().idxmin()
        f_found = float(freqs.loc[idx])

        rows.append({
            "measurement": str(sh).strip(),
            "freq_hz": f_found,
            "rho": float(df.loc[idx, SIP_RHO_COL]),
            "sigma_imag": float(df.loc[idx, SIP_IMAG_COL]),
            "freq_abs_err": float(abs(f_found - target_freq))
        })

    return rows


def is_run_root(p: Path) -> bool:
    return (p / ANALYSIS_DIR).is_dir() and (p / RAW_DIR).is_dir()


def process_one_run(run_root: Path, freq: float, sample_filter: str | None, dry_run: bool) -> bool:
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir = run_root / RAW_DIR
    results_dir = run_root / RESULTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")

    if sip_xlsx is None or log_xlsx is None:
        return False

    # Create Results folder only
    if not dry_run:
        results_dir.mkdir(parents=True, exist_ok=True)

    sat_map = read_logbook_saturation(log_xlsx, sample_filter=sample_filter)
    sip_rows = extract_sip_at_freq(sip_xlsx, target_freq=freq)

    # join
    joined = []
    for r in sip_rows:
        mid = r["measurement"]
        S = sat_map.get(mid, float("nan"))
        joined.append({
            "measurement": mid,
            "freq_hz": r["freq_hz"],
            "freq_abs_err": r["freq_abs_err"],
            "S": S,
            "rho": r["rho"],
            "sigma_imag": r["sigma_imag"]
        })

    out = pd.DataFrame(joined)
    out["logS"] = out["S"].apply(safe_log10)
    out["logRho"] = out["rho"].apply(safe_log10)
    out["logImag"] = out["sigma_imag"].apply(safe_log10)

    # fits on logs
    lin_r = fit_linear(out["logS"], out["logRho"])
    lin_i = fit_linear(out["logS"], out["logImag"])
    logi_r = fit_logistic(out["logS"], out["logRho"])
    logi_i = fit_logistic(out["logS"], out["logImag"])

    tag = freq_to_tag(freq)
    csv_name = f"S1_joined_{tag}.csv"
    png_name = f"S1_exponents_{tag}.png"
    txt_name = f"S1_fit_report_{tag}.txt"
    gp_name = "plot_exponents.gp"  # fixed name inside Results

    if dry_run:
        print(f"[DRY] {run_root}")
        print(f"  SIP: {sip_xlsx}")
        print(f"  LOG: {log_xlsx}")
        print(f"  Results: {results_dir}")
        return True

    # write CSV (stable column order)
    out2 = out[["measurement", "freq_hz", "freq_abs_err", "S", "rho", "sigma_imag", "logS", "logRho", "logImag"]]
    out2.to_csv(results_dir / csv_name, index=False)

    # write report
    def fmt(d, keys):
        return "\n".join([f"  {k}: {d.get(k)}" for k in keys])

    report = []
    report.append(f"Run root: {run_root}")
    report.append(f"SIP workbook: {sip_xlsx}")
    report.append(f"Logbook: {log_xlsx}")
    report.append(f"Target freq: {freq} Hz  (tag: {tag})")
    report.append("")
    report.append("Joined data (note: Degree Saturation converted to fraction if needed):")
    report.append(out2.to_string(index=False))
    report.append("")
    report.append("=== Linear fit on logs ===")
    report.append("Resistivity: log(rho) = a*log(S) + b   (Archie n = -a)")
    report.append(fmt(lin_r, ["n", "a", "b", "R2", "p_model", "n_obs"]) if False else
                 f"  a: {lin_r['a']:.6g}\n  b: {lin_r['b']:.6g}\n  n(=-a): {-lin_r['a']:.6g}\n  R2: {lin_r['R2']:.6g}\n  p_model: {lin_r['p_model']:.6g}\n  n_obs: {lin_r['n']}")
    report.append("")
    report.append("Imag cond: log(sig_imag) = a*log(S) + b   (saturation exponent p = a)")
    report.append(f"  a(p): {lin_i['a']:.6g}\n  b: {lin_i['b']:.6g}\n  R2: {lin_i['R2']:.6g}\n  p_model: {lin_i['p_model']:.6g}\n  n_obs: {lin_i['n']}")
    report.append("")
    report.append("=== Logistic fit on logs ===")
    report.append("Model: y = c + L / (1 + exp(-k*(x - x0)))   where x=log(S), y=log(property)")
    report.append("Resistivity logistic params:")
    report.append(f"  L: {logi_r['L']:.6g}\n  k: {logi_r['k']:.6g}\n  x0: {logi_r['x0']:.6g}\n  c: {logi_r['c']:.6g}\n"
                 f"  R2: {logi_r['R2']:.6g}\n  p_model: {logi_r['p_model']:.6g}\n"
                 f"  p_L: {logi_r['p_L']:.6g}\n  p_k: {logi_r['p_k']:.6g}\n  p_x0: {logi_r['p_x0']:.6g}\n  p_c: {logi_r['p_c']:.6g}\n"
                 f"  n_obs: {logi_r['n']}")
    report.append("")
    report.append("Imag cond logistic params:")
    report.append(f"  L: {logi_i['L']:.6g}\n  k: {logi_i['k']:.6g}\n  x0: {logi_i['x0']:.6g}\n  c: {logi_i['c']:.6g}\n"
                 f"  R2: {logi_i['R2']:.6g}\n  p_model: {logi_i['p_model']:.6g}\n"
                 f"  p_L: {logi_i['p_L']:.6g}\n  p_k: {logi_i['p_k']:.6g}\n  p_x0: {logi_i['p_x0']:.6g}\n  p_c: {logi_i['p_c']:.6g}\n"
                 f"  n_obs: {logi_i['n']}")

    (results_dir / txt_name).write_text("\n".join(report))

    # gnuplot + run
    write_gnuplot(results_dir, csv_name=csv_name, out_png=png_name, freq=freq,
                 lin_r=lin_r, logi_r=logi_r, lin_i=lin_i, logi_i=logi_i)

    # run gnuplot inside Results folder
    subprocess.run(["gnuplot", gp_name], cwd=str(results_dir), check=True)

    print(f"[OK] {run_root} -> {results_dir / png_name}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parent", help="Top folder that contains run folders (recursively scanned).")
    ap.add_argument("--freq", type=float, default=0.01, help="Target frequency in Hz (default: 0.01).")
    ap.add_argument("--sample", type=str, default=None,
                    help="Optional substring filter for Sample name in logbook (case-insensitive).")
    ap.add_argument("--dry-run", action="store_true", help="Do not create Results or run gnuplot.")
    args = ap.parse_args()

    base = Path(args.parent).expanduser().resolve()
    if not base.exists():
        raise SystemExit(f"Not found: {base}")

    count = 0
    for p in base.rglob("*"):
        if p.is_dir() and is_run_root(p):
            try:
                did = process_one_run(p, freq=args.freq, sample_filter=args.sample, dry_run=args.dry_run)
                if did:
                    count += 1
            except Exception as e:
                print(f"[FAIL] {p}: {e}")

    print(f"Done. Processed {count} run folder(s).")


if __name__ == "__main__":
    main()
