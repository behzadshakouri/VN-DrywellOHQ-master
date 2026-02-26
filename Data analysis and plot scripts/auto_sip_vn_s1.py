#!/usr/bin/env python3
import os
import sys
import math
import subprocess
from pathlib import Path

import pandas as pd

# -----------------------------
# User configuration (EDIT ME)
# -----------------------------

TARGET_FREQ_HZ = 0.01   # frequency to extract (Hz)
FREQ_TOL = 1e-9         # tolerance when matching frequency

# How to match "measurement id" (M1/M2/...) to a saturation value in the logbook
# You MUST set these based on your logbook structure.
LOGBOOK_SHEET = None  # None -> first sheet; or set to a sheet name like "S-1"
LOGBOOK_ID_COL = "Measurement"     # e.g. a column containing "M1", "M2", ...
LOGBOOK_SAT_COL = "Saturation"     # e.g. saturation as 0..1

# Output filenames
OUT_CSV = "S1_exponents_input_0p01Hz.csv"
GNUPLOT_SCRIPT = "plot_S1_exponents_0p01Hz.gp"
OUT_PNG = "VanNeys_S1_exponents_0p01Hz.png"

# Prefer these column names in the SIP output sheets (you can adapt if yours differ)
# If your first column is frequency, keep it.
SIP_FREQ_COL_CANDIDATES = ["Frequency", "Frequency (Hz)", "Freq", "Hz"]
SIP_RHO_COL_CANDIDATES  = ["Resistivity", "Resistivity (Ω·m)", "Resistivity (Ohm m)", "rho"]
SIP_IMAG_COL_CANDIDATES = ["Imaginary Conductivity", "Imaginary Conductivity (µS/cm)", "Imag", "sigma_imag"]

# File name patterns
ANALYSIS_XLSX_HINT = "output"     # in Analysis folder
LOGBOOK_XLSX_HINT  = "Log Book"   # in SIP_Raw_Data folder

# Folder names
ANALYSIS_DIR = "Analysis"
RAW_DIR = "SIP_Raw_Data"
PLOTS_DIR = "Plots"


# -----------------------------
# Helpers
# -----------------------------

def find_first_matching_file(folder: Path, must_contain: str, ext: str = ".xlsx") -> Path | None:
    if not folder.exists():
        return None
    for p in folder.rglob(f"*{ext}"):
        if must_contain.lower() in p.name.lower():
            return p
    return None

def pick_col(df: pd.DataFrame, candidates: list[str]) -> str:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    # fallback: return first numeric-ish column if candidates not found
    return ""

def to_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")

def safe_log10(x):
    x = float(x)
    if x <= 0 or math.isnan(x):
        return float("nan")
    return math.log10(x)

def linear_fit(x, y):
    """
    Returns slope, intercept, R^2 for y = slope*x + intercept
    """
    x = pd.Series(x).astype(float)
    y = pd.Series(y).astype(float)
    mask = (~x.isna()) & (~y.isna())
    x = x[mask]; y = y[mask]
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan")
    xmean = x.mean()
    ymean = y.mean()
    ss_xy = ((x - xmean) * (y - ymean)).sum()
    ss_xx = ((x - xmean) ** 2).sum()
    slope = ss_xy / ss_xx if ss_xx != 0 else float("nan")
    intercept = ymean - slope * xmean
    yhat = slope * x + intercept
    ss_res = ((y - yhat) ** 2).sum()
    ss_tot = ((y - ymean) ** 2).sum()
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else float("nan")
    return float(slope), float(intercept), float(r2)

def write_gnuplot(gp_path: Path, csv_name: str, out_png: str, fit_n: float, fit_r2: float, fit_p: float, fit_r2_p: float):
    # Two-panel plot like your image: resistivity and imag cond exponents
    gp = f"""\
reset
set datafile separator ','
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'

set multiplot layout 2,1

# --------------------
# Panel 1: Resistivity
# --------------------
set title 'Van Nuys S1: Saturation Exponent (Resistivity) at {TARGET_FREQ_HZ:g} Hz'
set grid
set key top right
set xlabel 'log(Saturation)'
set ylabel 'log(Resistivity Index)'
f1(x) = a1*x + b1
a1 = {-fit_n:.6g}
b1 = 0
plot '{csv_name}' using 5:6 with points pt 7 ps 1.6 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('{TARGET_FREQ_HZ:g}Hz: n = %.2f, R^2 = %.3f', -a1, {fit_r2:.6f})

# --------------------
# Panel 2: Imag cond
# --------------------
set title 'Van Nuys S1: Saturation Exponent (Imaginary Conductivity) at {TARGET_FREQ_HZ:g} Hz'
set grid
set key top right
set xlabel 'log(Saturation)'
set ylabel 'log(Imaginary Conductivity Index)'
f2(x) = a2*x + b2
a2 = {fit_p:.6g}
b2 = 0
plot '{csv_name}' using 5:7 with points pt 7 ps 1.6 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('{TARGET_FREQ_HZ:g}Hz: p = %.2f, R^2 = %.3f', a2, {fit_r2_p:.6f})

unset multiplot
"""
    gp_path.write_text(gp)


def process_one_run(run_root: Path) -> bool:
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir = run_root / RAW_DIR
    plots_dir = run_root / PLOTS_DIR

    sip_xlsx = find_first_matching_file(analysis_dir, ANALYSIS_XLSX_HINT, ".xlsx")
    log_xlsx = find_first_matching_file(raw_dir, LOGBOOK_XLSX_HINT, ".xlsx")

    if sip_xlsx is None or log_xlsx is None:
        return False

    plots_dir.mkdir(parents=True, exist_ok=True)

    # ---- Read logbook: build dict measurement_id -> saturation
    logbook = pd.read_excel(log_xlsx, sheet_name=LOGBOOK_SHEET)
    if LOGBOOK_ID_COL not in logbook.columns or LOGBOOK_SAT_COL not in logbook.columns:
        raise RuntimeError(
            f"[{run_root}] Logbook columns not found. "
            f"Expected '{LOGBOOK_ID_COL}' and '{LOGBOOK_SAT_COL}'. "
            f"Available columns: {list(logbook.columns)}"
        )
    sat_map = {}
    for _, r in logbook.iterrows():
        mid = str(r[LOGBOOK_ID_COL]).strip()
        s = to_float(r[LOGBOOK_SAT_COL])
        if mid and not math.isnan(s):
            sat_map[mid] = s

    # ---- Read SIP output: each sheet is a measurement (M1..)
    xls = pd.ExcelFile(sip_xlsx)
    rows = []
    for sh in xls.sheet_names:
        df = pd.read_excel(sip_xlsx, sheet_name=sh)

        # identify columns
        freq_col = pick_col(df, SIP_FREQ_COL_CANDIDATES)
        rho_col  = pick_col(df, SIP_RHO_COL_CANDIDATES)
        imag_col = pick_col(df, SIP_IMAG_COL_CANDIDATES)

        if not freq_col:
            # assume first column is frequency
            freq_col = df.columns[0]

        if not rho_col or not imag_col:
            # try positional fallbacks (based on your earlier structure)
            # columns: freq, imag_cond, phase, real_cond, resistivity
            if len(df.columns) >= 5:
                imag_col = df.columns[1]
                rho_col  = df.columns[4]
            else:
                continue

        # find the row closest to TARGET_FREQ_HZ
        freqs = df[freq_col].apply(to_float)
        idx = (freqs - TARGET_FREQ_HZ).abs().idxmin()
        if abs(freqs.loc[idx] - TARGET_FREQ_HZ) > max(FREQ_TOL, TARGET_FREQ_HZ*1e-6):
            # no close match
            continue

        rho = to_float(df.loc[idx, rho_col])
        sig_im = to_float(df.loc[idx, imag_col])

        # saturation from logbook map
        # default key: sheet name itself (e.g. "M1")
        mid = str(sh).strip()
        S = sat_map.get(mid, float("nan"))

        rows.append({
            "measurement": mid,
            "freq_hz": float(freqs.loc[idx]),
            "S": S,
            "rho": rho,
            "sigma_imag": sig_im,
        })

    if not rows:
        raise RuntimeError(f"[{run_root}] Found files but could not extract any rows at {TARGET_FREQ_HZ} Hz.")

    out = pd.DataFrame(rows)

    # compute logs used in your plots
    out["logS"] = out["S"].apply(safe_log10)
    out["logRho"] = out["rho"].apply(safe_log10)
    out["logImag"] = out["sigma_imag"].apply(safe_log10)

    # fit log(rho) vs log(S): slope = -n
    slope_r, intercept_r, r2_r = linear_fit(out["logS"], out["logRho"])
    n = -slope_r

    # fit log(sig_im) vs log(S): slope = p
    slope_i, intercept_i, r2_i = linear_fit(out["logS"], out["logImag"])
    p = slope_i

    # write CSV into Plots
    csv_path = plots_dir / OUT_CSV
    # columns in order so gnuplot can use fixed indices if you want
    out2 = out[["measurement","freq_hz","S","rho","logS","logRho","logImag"]]
    out2.to_csv(csv_path, index=False)

    # write gnuplot
    gp_path = plots_dir / GNUPLOT_SCRIPT
    write_gnuplot(
        gp_path=gp_path,
        csv_name=OUT_CSV,
        out_png=OUT_PNG,
        fit_n=n,
        fit_r2=r2_r,
        fit_p=p,
        fit_r2_p=r2_i
    )

    # run gnuplot inside Plots folder (so relative paths work)
    subprocess.run(["gnuplot", gp_path.name], cwd=str(plots_dir), check=True)
    print(f"[OK] {run_root} -> {plots_dir/OUT_PNG}")
    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 auto_sip_vn_s1.py /path/to/parent_folder")
        sys.exit(2)

    base = Path(sys.argv[1]).expanduser().resolve()
    if not base.exists():
        print(f"Not found: {base}")
        sys.exit(2)

    count = 0
    # A "run root" is any folder containing Analysis/ and SIP_Raw_Data/
    for p in base.rglob("*"):
        if p.is_dir() and (p/ANALYSIS_DIR).is_dir() and (p/RAW_DIR).is_dir():
            try:
                did = process_one_run(p)
                if did:
                    count += 1
            except Exception as e:
                print(f"[FAIL] {p}: {e}")

    print(f"Done. Processed {count} run folder(s).")

if __name__ == "__main__":
    main()
