#!/usr/bin/env python3
"""
Fully automatic SIP + saturation pipeline:

For each "run root" directory that contains:
  - Analysis/      (contains *_output*.xlsx with sheets M1, M2, ...)
  - SIP_Raw_Data/  (contains *Log Book*.xlsx with a "Saturation" sheet)
  - Plots/         (will be created if missing)

This script will:
  1) Read degree of saturation S for each measurement (M1, M2, ...) from the logbook
  2) Read resistivity rho and imaginary conductivity sigma_imag at a target frequency from the SIP output workbook
  3) Join them on measurement ID
  4) Write a CSV into Plots/
  5) Generate and run a gnuplot script to make a 2-panel PNG like your figure

Usage:
  python3 auto_sip_pipeline.py /path/to/parent --freq 0.01

Notes:
  - The script is robust to logbooks where 'Sample' is only filled in the first row:
    it forward-fills 'Sample' so every row has a sample label.
  - If multiple samples exist in one logbook, it automatically chooses the sample group
    with the best overlap with the measurement sheet names (M1..).
"""
import argparse
import math
import subprocess
from pathlib import Path

import pandas as pd


# -----------------------------
# Defaults (match your folders)
# -----------------------------
ANALYSIS_DIR = "Analysis"
RAW_DIR = "SIP_Raw_Data"
PLOTS_DIR = "Plots"

ANALYSIS_XLSX_HINT = "output"      # file name contains this (case-insensitive)
LOGBOOK_XLSX_HINT  = "Log Book"    # file name contains this (case-insensitive)

LOGBOOK_SHEET = "Saturation"
LOGBOOK_SAMPLE_COL = "Sample"
LOGBOOK_ID_COL = "Measurements"
LOGBOOK_SAT_COL = "Degree Saturation M2"   # in percent (0..100) in your uploaded file

# SIP output columns (exact headers in your uploaded Van Nuys S1_output.xlsx)
SIP_FREQ_COL = "Frequency (Hz)"
SIP_RHO_COL  = "Resistivity (Ω·m)"
SIP_IMAG_COL = "Imaginary_Conductivity (µS/cm)"


def find_first_excel(folder: Path, name_hint: str) -> Path | None:
    if not folder.exists():
        return None
    hint = name_hint.lower()
    for p in folder.rglob("*.xlsx"):
        if hint in p.name.lower():
            return p
    return None


def safe_log10(x: float) -> float:
    try:
        x = float(x)
    except Exception:
        return float("nan")
    if not math.isfinite(x) or x <= 0:
        return float("nan")
    return math.log10(x)


def linear_fit(x, y):
    """
    y = m*x + b ; return m, b, R^2 using ordinary least squares
    """
    x = pd.Series(x, dtype="float64")
    y = pd.Series(y, dtype="float64")
    mask = (~x.isna()) & (~y.isna())
    x = x[mask]; y = y[mask]
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan")
    xm = x.mean(); ym = y.mean()
    ss_xx = ((x - xm) ** 2).sum()
    ss_xy = ((x - xm) * (y - ym)).sum()
    m = ss_xy / ss_xx if ss_xx != 0 else float("nan")
    b = ym - m * xm
    yhat = m * x + b
    ss_res = ((y - yhat) ** 2).sum()
    ss_tot = ((y - ym) ** 2).sum()
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else float("nan")
    return float(m), float(b), float(r2)


def choose_best_sample_group(log_df: pd.DataFrame, sheet_ids: set[str]) -> pd.DataFrame:
    """
    If multiple samples exist, pick the sample whose Measurements overlap most with sheet_ids.
    """
    if LOGBOOK_SAMPLE_COL not in log_df.columns:
        return log_df

    tmp = log_df.copy()
    tmp[LOGBOOK_SAMPLE_COL] = tmp[LOGBOOK_SAMPLE_COL].ffill()

    if tmp[LOGBOOK_SAMPLE_COL].isna().all():
        return log_df

    best = None
    best_score = -1
    for sample, g in tmp.groupby(LOGBOOK_SAMPLE_COL, dropna=True):
        mids = set(str(v).strip() for v in g.get(LOGBOOK_ID_COL, pd.Series([], dtype=object)).dropna().tolist())
        score = len(mids & sheet_ids)
        if score > best_score:
            best_score = score
            best = g

    return best if best is not None else tmp


def write_gnuplot(gp_path: Path, csv_name: str, out_png: str, freq: float, n: float, r2_n: float, p: float, r2_p: float):
    # We use columns by name via `using` with `column("name")` for safety.
    gp = f"""\
reset
set datafile separator ','
set term pngcairo size 1200,1600 enhanced font 'Arial,26'
set output '{out_png}'

set multiplot layout 2,1

# --------------------
# Panel 1: Resistivity
# --------------------
set title sprintf('Van Nuys S1: Saturation Exponent (Resistivity) at %.5g Hz', {freq})
set grid
set key top right
set xlabel 'log(Saturation)'
set ylabel 'log(Resistivity Index)'

f1(x) = a1*x + b1
a1 = {-n:.12g}     # because log(rho) = -n*log(S) + C  => slope = -n
b1 = 0

plot '{csv_name}' using (column("logS")):(column("logRho")) with points pt 7 ps 1.6 title 'data', \\
     f1(x) with lines dt 2 lw 3 title sprintf('%.5gHz: n = %.2f, R^2 = %.3f', {freq}, {n:.6g}, {r2_n:.6g})

# --------------------
# Panel 2: Imag cond
# --------------------
set title sprintf('Van Nuys S1: Saturation Exponent (Imaginary Conductivity) at %.5g Hz', {freq})
set grid
set key top right
set xlabel 'log(Saturation)'
set ylabel 'log(Imaginary Conductivity Index)'

f2(x) = a2*x + b2
a2 = {p:.12g}
b2 = 0

plot '{csv_name}' using (column("logS")):(column("logImag")) with points pt 7 ps 1.6 title 'data', \\
     f2(x) with lines dt 2 lw 3 title sprintf('%.5gHz: p = %.2f, R^2 = %.3f', {freq}, {p:.6g}, {r2_p:.6g})

unset multiplot
"""
    gp_path.write_text(gp)


def process_run_root(run_root: Path, freq: float, force_sample: str | None, dry_run: bool) -> bool:
    analysis_dir = run_root / ANALYSIS_DIR
    raw_dir = run_root / RAW_DIR
    plots_dir = run_root / PLOTS_DIR

    sip_xlsx = find_first_excel(analysis_dir, ANALYSIS_XLSX_HINT)
    log_xlsx = find_first_excel(raw_dir, LOGBOOK_XLSX_HINT)

    if sip_xlsx is None or log_xlsx is None:
        return False

    # Load SIP workbook
    sip_book = pd.ExcelFile(sip_xlsx)
    sheet_ids = set(sip_book.sheet_names)

    # Load logbook saturation sheet
    log_df = pd.read_excel(log_xlsx, sheet_name=LOGBOOK_SHEET)

    # forward-fill sample, and optionally filter
    if LOGBOOK_SAMPLE_COL in log_df.columns:
        log_df[LOGBOOK_SAMPLE_COL] = log_df[LOGBOOK_SAMPLE_COL].ffill()

    if force_sample:
        # keep rows whose sample string contains the force_sample token (case-insensitive)
        if LOGBOOK_SAMPLE_COL in log_df.columns:
            mask = log_df[LOGBOOK_SAMPLE_COL].astype(str).str.lower().str.contains(force_sample.lower(), na=False)
            if mask.any():
                log_df = log_df[mask]
        # else: ignore (no sample column)

    # If multiple samples remain, choose best overlap
    log_df = choose_best_sample_group(log_df, sheet_ids)

    # Build saturation mapping (Measurements -> S in 0..1)
    if (LOGBOOK_ID_COL not in log_df.columns) or (LOGBOOK_SAT_COL not in log_df.columns):
        raise RuntimeError(
            f"[{run_root}] Logbook missing expected columns.\n"
            f"Need: '{LOGBOOK_ID_COL}' and '{LOGBOOK_SAT_COL}' in sheet '{LOGBOOK_SHEET}'.\n"
            f"Have: {list(log_df.columns)}"
        )

    sat_map = {}
    for _, r in log_df.iterrows():
        mid = str(r[LOGBOOK_ID_COL]).strip()
        s = r[LOGBOOK_SAT_COL]
        try:
            s = float(s)
        except Exception:
            s = float("nan")
        if not mid or not math.isfinite(s):
            continue
        # In your logbook this appears as percent; convert if needed
        if s > 1.5:
            s = s / 100.0
        sat_map[mid] = s

    # Extract per-sheet values at target frequency
    rows = []
    for sh in sip_book.sheet_names:
        df = pd.read_excel(sip_xlsx, sheet_name=sh)

        # Basic column checks (fall back to positional if headers differ)
        if SIP_FREQ_COL not in df.columns:
            freq_col = df.columns[0]
        else:
            freq_col = SIP_FREQ_COL

        rho_col = SIP_RHO_COL if SIP_RHO_COL in df.columns else (df.columns[4] if len(df.columns) >= 5 else None)
        imag_col = SIP_IMAG_COL if SIP_IMAG_COL in df.columns else (df.columns[1] if len(df.columns) >= 2 else None)
        if rho_col is None or imag_col is None:
            continue

        freqs = pd.to_numeric(df[freq_col], errors="coerce")
        if freqs.isna().all():
            continue

        idx = (freqs - freq).abs().idxmin()
        fval = float(freqs.loc[idx])

        rho = pd.to_numeric(df.loc[idx, rho_col], errors="coerce")
        imag = pd.to_numeric(df.loc[idx, imag_col], errors="coerce")

        S = sat_map.get(str(sh).strip(), float("nan"))
        rows.append({
            "Measurement": str(sh).strip(),
            "freq_hz": fval,
            "S": float(S) if math.isfinite(S) else float("nan"),
            "rho": float(rho) if math.isfinite(rho) else float("nan"),
            "sigma_imag": float(imag) if math.isfinite(imag) else float("nan"),
        })

    out = pd.DataFrame(rows)
    # keep only rows with saturation and positive electrical values
    out = out[(out["S"] > 0) & (out["rho"] > 0) & (out["sigma_imag"] > 0)].copy()
    if out.empty:
        raise RuntimeError(f"[{run_root}] No valid joined rows. Check measurement IDs and saturation mapping.")

    # logs
    out["logS"] = out["S"].apply(safe_log10)
    out["logRho"] = out["rho"].apply(safe_log10)
    out["logImag"] = out["sigma_imag"].apply(safe_log10)

    # fits
    m_r, b_r, r2_r = linear_fit(out["logS"], out["logRho"])
    n = -m_r
    m_i, b_i, r2_i = linear_fit(out["logS"], out["logImag"])
    p = m_i

    # outputs
    plots_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{freq:g}".replace(".", "p")
    csv_name = f"S1_exponents_input_{tag}Hz.csv"
    gp_name = f"plot_S1_exponents_{tag}Hz.gp"
    png_name = f"Van_Nuys_S1_exponents_{tag}Hz.png"

    csv_path = plots_dir / csv_name
    gp_path = plots_dir / gp_name

    out.to_csv(csv_path, index=False)

    write_gnuplot(gp_path, csv_name, png_name, freq, n, r2_r, p, r2_i)

    if dry_run:
        print(f"[DRY] {run_root}\n  sip: {sip_xlsx}\n  log: {log_xlsx}\n  -> {csv_path.name}, {gp_path.name}, {png_name}")
        return True

    # Run gnuplot in Plots dir so relative paths work
    subprocess.run(["gnuplot", gp_path.name], cwd=str(plots_dir), check=True)
    print(f"[OK] {run_root} -> {plots_dir / png_name}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base_dir", help="Parent directory to scan recursively")
    ap.add_argument("--freq", type=float, default=0.01, help="Target frequency in Hz (default 0.01)")
    ap.add_argument("--sample", type=str, default=None,
                    help="Optional token to select Sample rows in logbook (e.g. 'Van Nuys S-1'). Case-insensitive substring match.")
    ap.add_argument("--dry-run", action="store_true", help="Do not run gnuplot; just write CSV/GP and report.")
    args = ap.parse_args()

    base = Path(args.base_dir).expanduser().resolve()
    if not base.exists():
        raise SystemExit(f"Not found: {base}")

    processed = 0
    for p in base.rglob("*"):
        if p.is_dir() and (p / ANALYSIS_DIR).is_dir() and (p / RAW_DIR).is_dir():
            try:
                did = process_run_root(p, args.freq, args.sample, args.dry_run)
                if did:
                    processed += 1
            except Exception as e:
                print(f"[FAIL] {p}: {e}")

    print(f"Done. Processed {processed} run folder(s).")


if __name__ == "__main__":
    main()
