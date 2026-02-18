import argparse
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

EXCEL_EPOCH = datetime(1899, 12, 30)

def excel_serial_to_datetime(x: float) -> datetime:
    return EXCEL_EPOCH + timedelta(days=float(x))

def datetime_to_excel_serial(dt) -> float:
    dt = pd.to_datetime(dt).to_pydatetime()
    return (dt - EXCEL_EPOCH).total_seconds() / 86400.0

def as_time_both(v):
    """
    Return (time_excel, time_str) from a cell value v.
    Handles:
      - numeric Excel serial (e.g., 45763.588194)
      - pandas Timestamp / datetime
      - string dates
    """
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return (np.nan, "UNKNOWN")

    # numeric Excel serial directly
    if isinstance(v, (int, float, np.integer, np.floating)):
        x = float(v)
        if 20000.0 <= x <= 80000.0:  # sanity range
            dt = excel_serial_to_datetime(x)
            return (x, dt.strftime("%Y-%m-%d %H:%M"))

    # parse datetime/string
    try:
        t = pd.to_datetime(v, errors="coerce")
        if pd.isna(t):
            return (np.nan, "UNKNOWN")
        x = datetime_to_excel_serial(t)
        return (x, t.strftime("%Y-%m-%d %H:%M"))
    except Exception:
        return (np.nan, "UNKNOWN")

def nearest_time_cell_value(row0, c):
    """Find nearest non-empty timestamp cell to the LEFT of column c in row0."""
    for cc in range(c, -1, -1):
        v = row0.iloc[cc]
        if pd.notna(v):
            return v
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xlsx", help="Input .xlsx file")
    ap.add_argument("--sheet", default=None, help="Sheet name (default: first sheet)")
    ap.add_argument("--borehole", default="ERT-3", help="Borehole label to store in output")
    ap.add_argument("--out", default="ert_tidy.csv", help="Output CSV")
    ap.add_argument("--depth_header", default="Depth m", help="Header text for depth column")
    ap.add_argument("--sm_header", default="Soil Moisture %", help="Header text for soil moisture column")
    ap.add_argument("--header_row", type=int, default=1, help="0-based row index containing column headers")
    ap.add_argument("--time_row", type=int, default=0, help="0-based row index containing time labels")
    args = ap.parse_args()

    xl = pd.ExcelFile(args.xlsx)
    sheet = args.sheet or xl.sheet_names[0]
    raw = pd.read_excel(args.xlsx, sheet_name=sheet, header=None)

    if raw.shape[0] <= max(args.header_row, args.time_row):
        raise SystemExit(f"Sheet too short: shape={raw.shape}. Check --header_row/--time_row.")

    hdr = raw.iloc[args.header_row]
    row0 = raw.iloc[args.time_row]

    # Find header positions
    depth_cols = [c for c in range(raw.shape[1]) if str(hdr.iloc[c]).strip() == args.depth_header]
    sm_cols    = [c for c in range(raw.shape[1]) if str(hdr.iloc[c]).strip() == args.sm_header]

    if not depth_cols or not sm_cols:
        raise SystemExit(
            f"Could not find headers. Found Depth cols={depth_cols}, SM cols={sm_cols}. "
            f"Check --depth_header/--sm_header/--header_row."
        )

    # Pair depth + moisture columns by nearest moisture column to the right
    pairs = []
    for dc in depth_cols:
        candidates = [sc for sc in sm_cols if sc >= dc]
        if candidates:
            sc = min(candidates, key=lambda x: abs(x - dc))
            pairs.append((dc, sc))

    out_rows = []
    data_start = args.header_row + 1

    for (dc, sc) in pairs:
        time_cell = nearest_time_cell_value(row0, dc)
        t_excel, t_str = as_time_both(time_cell)

        rr = data_start
        while rr < raw.shape[0]:
            d = raw.iat[rr, dc]
            m = raw.iat[rr, sc]

            # stop at first fully-empty row for this block
            if pd.isna(d) and pd.isna(m):
                break

            if pd.notna(d) and pd.notna(m):
                try:
                    out_rows.append({
                        "borehole": args.borehole,
                        "time_str": t_str,
                        "time_excel": float(t_excel) if pd.notna(t_excel) else np.nan,
                        "depth_m": float(d),
                        "soil_moisture_pct": float(m),
                    })
                except Exception:
                    pass

            rr += 1

    tidy = pd.DataFrame(out_rows)
    if tidy.empty:
        raise SystemExit(
            "No data extracted (empty output).\n"
            "Likely causes: wrong --header_row/--time_row, merged cells, or headers differ.\n"
        )

    # Force column order (guaranteed)
    tidy = tidy[["borehole", "time_str", "time_excel", "depth_m", "soil_moisture_pct"]]

    tidy.to_csv(args.out, index=False)
    print(f"Wrote {args.out}  (rows={len(tidy)})")

if __name__ == "__main__":
    main()
