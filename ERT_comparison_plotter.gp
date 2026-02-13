#!/usr/bin/gnuplot
# ============================================================
# ERT vs Model snapshots (single run, CURRENT DIRECTORY)
# Works with BOTH filename conventions:
#   NEW: ERTsnap-t<ttok>-<bh>.csv
#   OLD: ERTsnap_t<ttok>_<bh>.csv
#
# CSV columns:
#   depth_m,theta_model,theta_obs,time
#
# Outputs:
# (A) One figure per (borehole, time): model line + obs points
# (B) One figure per borehole: 5 timesteps overlaid
# (C) One figure per borehole: 5 subplots (1x5)
# ============================================================

# -------------------------------
# Global appearance
# -------------------------------
set datafile separator ','
set grid
set tics out
set border lw 2

# Depth increases downward
set ylabel "Depth (m)" font "Arial,32"
set yrange [*:*] reverse

# If model is 0..1 and obs is 0..100, set MODEL_SCALE=100.
# If both are same units, set MODEL_SCALE=1.
MODEL_SCALE = 100

# -------------------------------
# Styles
# -------------------------------
# Model line (black)
set style line 11 lc rgb "#000000" lw 3.5 dt 1
# Obs points (red)
set style line 12 lc rgb "#d60000" lw 2.6 pt 7 ps 1.15

# For multi-time overlay: 5 model dash styles (black)
set style line 21 lc rgb "#000000" lw 3.5 dt 1
set style line 22 lc rgb "#000000" lw 3.5 dt 2
set style line 23 lc rgb "#000000" lw 3.5 dt 3
set style line 24 lc rgb "#000000" lw 3.5 dt 4
set style line 25 lc rgb "#000000" lw 3.5 dt 5

# For multi-time overlay: 5 obs point styles (red)
set style line 31 lc rgb "#d60000" lw 2.4 pt 7 ps 1.10 dt 1
set style line 32 lc rgb "#d60000" lw 2.4 pt 7 ps 1.10 dt 2
set style line 33 lc rgb "#d60000" lw 2.4 pt 7 ps 1.10 dt 3
set style line 34 lc rgb "#d60000" lw 2.4 pt 7 ps 1.10 dt 4
set style line 35 lc rgb "#d60000" lw 2.4 pt 7 ps 1.10 dt 5

# -------------------------------
# Inputs
# -------------------------------
BHs  = "ERT-3 ERT-5"

TTOK = "45763p588194 45763p682639 45763p770833 45764p673611 45764p729167"

# IMPORTANT: avoid '_' here (gnuplot enhanced treats _ as subscript)
TLBL = "2025-04-16-14-07 2025-04-16-16-23 2025-04-16-18-30 2025-04-17-16-10 2025-04-17-17-30"

# NEW naming pieces
PREFIX_NEW = "ERTsnap-t"
MID_NEW    = "-"
# OLD naming pieces
PREFIX_OLD = "ERTsnap_t"
MID_OLD    = "_"
SUFFIX     = ".csv"

# -------------------------------
# Helpers
# -------------------------------
stripnl(s) = (strlen(s) > 0 && s[strlen(s):strlen(s)] eq "\n") ? s[1:strlen(s)-1] : s

# Find file in current directory, try NEW then OLD
# NEW: ERTsnap-t<ttok>-<bh>.csv
# OLD: ERTsnap_t<ttok>_<bh>.csv
find_file(bh, ttok) = stripnl(system(sprintf( \
  "bash -lc 'f1=\"%s%s%s%s%s\"; f2=\"%s%s%s%s%s\"; " . \
  "if [ -f \"$f1\" ]; then echo \"$f1\"; " . \
  "elif [ -f \"$f2\" ]; then echo \"$f2\"; " . \
  "else echo \"\"; fi'", \
  PREFIX_NEW, ttok, MID_NEW, bh, SUFFIX, \
  PREFIX_OLD, ttok, MID_OLD, bh, SUFFIX)))

# When a file is missing, draw a harmless dummy point + label (no NaN crash)
plot_missing(msg) = sprintf( \
  "unset label; set label 1 '%s' at graph 0.5, graph 0.5 center font 'Arial,28'; " . \
  "plot '-' using 1:2 with points pt 7 ps 0.01 notitle; 0 0; e", msg)

# ============================================================
# (A) One plot per (borehole,time): model line + obs points
# ============================================================
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set key right top
set xlabel "Soil Moisture" font "Arial,32"

do for [ib=1:words(BHs)] {
  bh = word(BHs, ib)

  do for [it=1:words(TTOK)] {
    ttok = word(TTOK, it)
    tlbl = word(TLBL, it)

    f = find_file(bh, ttok)

    outfile = sprintf("ERT_cmp_%s_%s.png", bh, tlbl)
    set output outfile
    set title sprintf("%s | %s | Model (line) vs Obs (points)", bh, tlbl) font "Arial,32"

    if (strlen(f) == 0) {
      print sprintf("WARNING: Missing snapshot file for (bh=%s ttok=%s)", bh, ttok)
      eval(plot_missing(sprintf("MISSING\\n%s\\n%s", bh, ttok)))
    } else {
      unset label
      plot \
        f every ::1 using (($2)*MODEL_SCALE):1 with lines  ls 11 title "Model", \
        f every ::1 using ($3):1               with points ls 12 title "Obs"
    }

    unset output
  }
}

# ============================================================
# (B) Per borehole: 5 timesteps overlaid (model lines + obs points)
# ============================================================
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set key right top
set xlabel "Soil Moisture" font "Arial,32"

do for [ib=1:words(BHs)] {
  bh = word(BHs, ib)

  set output sprintf("ERT_profile_%s_all5_overlaid.png", bh)
  set title sprintf("%s | 5 ERT times (Model lines, Obs points)", bh) font "Arial,32"

  plotcmd = ""
  sep = ""

  do for [it=1:words(TTOK)] {
    ttok = word(TTOK, it)
    tlbl = word(TLBL, it)
    f = find_file(bh, ttok)

    if (strlen(f) == 0) {
      print sprintf("WARNING: Missing snapshot for overlay: bh=%s ttok=%s", bh, ttok)
    } else {
      plotcmd = plotcmd . sep . \
        sprintf("'%s' every ::1 using (($2)*MODEL_SCALE):1 with lines  ls %d title 'Model %s'", f, 20+it, tlbl) . ", " . \
        sprintf("'%s' every ::1 using ($3):1               with points ls %d title 'Obs  %s'",  f, 30+it, tlbl)
      sep = ", "
    }
  }

  if (strlen(plotcmd) > 0) {
    unset label
    eval("plot " . plotcmd)
  } else {
    eval(plot_missing(sprintf("NO FILES\\n%s", bh)))
  }

  unset output
}

# ============================================================
# (C) Per borehole: 5 subplots side-by-side (1x5)
# ============================================================
set terminal pngcairo enhanced font "Arial,22" size 2600,650
unset key

do for [ib=1:words(BHs)] {
  bh = word(BHs, ib)

  set output sprintf("ERT_profile_%s_5panels.png", bh)
  set multiplot layout 1,5 rowsfirst

  do for [it=1:words(TTOK)] {
    ttok = word(TTOK, it)
    tlbl = word(TLBL, it)
    f = find_file(bh, ttok)

    # y labels only on first panel
    if (it == 1) {
      set ytics
      set ylabel "Depth (m)" font "Arial,26"
    } else {
      unset ytics
      set ylabel ""
    }

    set xlabel "Moisture" font "Arial,24"
    set title tlbl font "Arial,24"

    if (strlen(f) == 0) {
      print sprintf("WARNING: Missing snapshot for panels: bh=%s ttok=%s", bh, ttok)
      eval(plot_missing(sprintf("MISSING\\n%s\\n%s", bh, ttok)))
    } else {
      unset label
      plot \
        f every ::1 using (($2)*MODEL_SCALE):1 with lines  lc rgb "#000000" lw 3.0 dt 1 notitle, \
        f every ::1 using ($3):1               with points lc rgb "#d60000" lw 2.2 pt 7 ps 1.0 notitle
    }
  }

  unset multiplot
  unset output
}

unset output

