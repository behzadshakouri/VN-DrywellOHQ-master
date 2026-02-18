#!/usr/bin/gnuplot
# ============================================================
# ERT vs Model snapshots (single run, CURRENT DIRECTORY)
# Works with BOTH filename conventions:
#   NEW: ERTsnap-t<ttok>-<bh>.csv
#   OLD: ERTsnap_t<ttok>_<bh>.csv
#
# Snapshot CSV columns:
#   depth_m,theta_model,theta_obs,time
#
# OBS source SWITCH:
#   USE_TIDY_OBS = 1  -> read obs from <bh>_tidy.csv
#                       columns: borehole,time_str,time_excel,depth_m,soil_moisture_pct
#   USE_TIDY_OBS = 0  -> old behavior: obs from snapshot col 3
#
# Outputs:
# (A) One figure per (borehole, time): model line+points + obs line+points
# (B) One figure per borehole: 5 timesteps overlaid
# (C) One figure per borehole: 5 subplots (1x5)
# ============================================================

# -------------------------------
# SWITCH
# -------------------------------
USE_TIDY_OBS = 1

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
# Model (black line + points; different marker shape)
set style line 11 lc rgb "#000000" lw 3.5 dt 1 pt 5 ps 1.10

# Observed (red line + points)
set style line 12 lc rgb "#d60000" lw 3.0 dt 1 pt 7 ps 1.20

# Overlay: 5 model line+point dash styles (black)
set style line 21 lc rgb "#000000" lw 4 dt 1 pt 5 ps 1.10
set style line 22 lc rgb "#000000" lw 4 dt 2 pt 5 ps 1.10
set style line 23 lc rgb "#000000" lw 4 dt 3 pt 5 ps 1.10
set style line 24 lc rgb "#000000" lw 4 dt 4 pt 5 ps 1.10
set style line 25 lc rgb "#000000" lw 4 dt 5 pt 5 ps 1.10

# Overlay: 5 obs line+point styles (red)
set style line 31 lc rgb "#d60000" lw 3.0 dt 1 pt 7 ps 1.10
set style line 32 lc rgb "#d60000" lw 3.0 dt 2 pt 7 ps 1.10
set style line 33 lc rgb "#d60000" lw 3.0 dt 3 pt 7 ps 1.10
set style line 34 lc rgb "#d60000" lw 3.0 dt 4 pt 7 ps 1.10
set style line 35 lc rgb "#d60000" lw 3.0 dt 5 pt 7 ps 1.10

# -------------------------------
# Inputs
# -------------------------------
BHs  = "ERT-3 ERT-5"
TTOK = "45763p588194 45763p682639 45763p770833 45764p673611 45764p729167"

# For DISPLAY (legend/title): use underscore + colon, then format to "YYYY-MM-DD HH:MM"
# (no gsub() in gnuplot; we format by slicing)
TLBL_RAW  = "2025-04-16_14:07 2025-04-16_16:23 2025-04-16_18:30 2025-04-17_16:10 2025-04-17_17:30"

# For FILENAMES (safe; no colon)
TLBL_FILE = "2025-04-16-14-07 2025-04-16-16-23 2025-04-16-18-30 2025-04-17-16-10 2025-04-17-17-30"

# Snapshot naming pieces
PREFIX_NEW = "ERTsnap-t"
MID_NEW    = "-"
PREFIX_OLD = "ERTsnap_t"
MID_OLD    = "_"
SUFFIX     = ".csv"

# -------------------------------
# Helpers
# -------------------------------
stripnl(s) = (strlen(s) > 0 && s[strlen(s):strlen(s)] eq "\n") ? s[1:strlen(s)-1] : s

# Turn "YYYY-MM-DD_HH:MM" into "YYYY-MM-DD HH:MM" without gsub()
fmt_time(s) = s[1:10] . " " . s[12:16]

# Find snapshot file in current directory, try NEW then OLD
find_file(bh, ttok) = stripnl(system(sprintf( \
  "sh -c 'f1=\"%s%s%s%s%s\"; f2=\"%s%s%s%s%s\"; " . \
  "if [ -f \"$f1\" ]; then echo \"$f1\"; " . \
  "elif [ -f \"$f2\" ]; then echo \"$f2\"; " . \
  "else echo \"\"; fi'", \
  PREFIX_NEW, ttok, MID_NEW, bh, SUFFIX, \
  PREFIX_OLD, ttok, MID_OLD, bh, SUFFIX)))

# When a file is missing, draw a harmless dummy point + label (no NaN crash)
plot_missing(msg) = sprintf( \
  "unset label; set label 1 '%s' at graph 0.5, graph 0.5 center font 'Arial,28'; " . \
  "plot '-' using 1:2 with points pt 7 ps 0.01 notitle; 0 0; e", msg)

# -------------------------------
# Tidy OBS helpers (temp file, comma-separated)
#   tidy file: <bh>_tidy.csv  in CURRENT DIRECTORY
#   columns: borehole,time_str,time_excel,depth_m,soil_moisture_pct
# -------------------------------
tidy_file(bh) = sprintf("%s_tidy.csv", bh)
obs_tmp_name(bh, ttok) = sprintf("__obs_%s_%s.dat", bh, ttok)

# tol=1e-5 to match rounded TTOK tokens vs Excel repeating decimals
make_obs_tmp(bh, ttok) = sprintf( \
  "sh -c 'in=\"%s\"; out=\"%s\"; " . \
  "if [ ! -f \"$in\" ]; then : > \"$out\"; exit 0; fi; " . \
  "awk -F, -v tok=\"%s\" " . \
  "\"BEGIN{gsub(/p/,\\\".\\\",tok); t=tok+0; tol=1e-5} " . \
  "NR==1{next} " . \
  "{x=\\$3; gsub(/\\\"/,\\\"\\\",x); x+=0; " . \
  " m=\\$5; d=\\$4; gsub(/\\\"/,\\\"\\\",m); gsub(/\\\"/,\\\"\\\",d); " . \
  " m+=0; d+=0; " . \
  " if (x>=t-tol && x<=t+tol) { if (m==m && d==d) print m \\\",\\\" d; } }\" " . \
  "\"$in\" > \"$out\"' ", \
  tidy_file(bh), obs_tmp_name(bh,ttok), ttok)

obs_tmp_count(bh, ttok) = int(system(sprintf( \
  "sh -c 'f=\"%s\"; if [ -f \"$f\" ]; then wc -l < \"$f\"; else echo 0; fi'", \
  obs_tmp_name(bh,ttok))))

# ============================================================
# (A) One plot per (borehole,time): model line+points + obs line+points
# ============================================================
set terminal pngcairo enhanced font "Arial,28" size 1200,800
set key right top
set xlabel "Soil Moisture" font "Arial,32"

do for [ib=1:words(BHs)] {
  bh = word(BHs, ib)

  do for [it=1:words(TTOK)] {
    ttok = word(TTOK, it)
    tlbl_disp = fmt_time(word(TLBL_RAW, it))
    tlbl_file = word(TLBL_FILE, it)

    f = find_file(bh, ttok)

    outfile = sprintf("ERT_cmp_%s_%s.png", bh, tlbl_file)
    set output outfile
    set title sprintf("%s | %s | Model vs Obs", bh, tlbl_disp) font "Arial,32"

    if (strlen(f) == 0) {
      print sprintf("WARNING: Missing snapshot file for (bh=%s ttok=%s)", bh, ttok)
      eval(plot_missing(sprintf("MISSING\\n%s\\n%s", bh, ttok)))
    } else {
      unset label

      if (USE_TIDY_OBS) {
        system(make_obs_tmp(bh, ttok))
        nobs = obs_tmp_count(bh, ttok)

        if (nobs <= 0) {
          print sprintf("WARNING: No tidy OBS rows (bh=%s ttok=%s) in %s", bh, ttok, tidy_file(bh))
          plot \
            f every ::1 using (($2)*MODEL_SCALE):1 with linespoints ls 11 title "Model"
        } else {
          plot \
            f every ::1 using (($2)*MODEL_SCALE):1 with linespoints ls 11 title "Model", \
            obs_tmp_name(bh, ttok) using 1:2        with linespoints ls 12 title "Obs"
        }
      } else {
        plot \
          f every ::1 using (($2)*MODEL_SCALE):1 with linespoints ls 11 title "Model", \
          f every ::1 using ($3):1               with linespoints ls 12 title "Obs"
      }
    }

    unset output
  }
}

# ============================================================
# (B) Per borehole: 5 timesteps overlaid
# ============================================================
set terminal pngcairo enhanced font "Arial,28" size 1400,800
set xlabel "Soil Moisture" font "Arial,32"

do for [ib=1:words(BHs)] {
  bh = word(BHs, ib)

  set output sprintf("ERT_profile_%s_all5_overlaid.png", bh)
  set title sprintf("%s | 5 ERT times (Model vs Obs)", bh) font "Arial,32"

  set key outside right top
  set key font "Arial,22"
  set key samplen 2 spacing 1.15
  set key opaque box

  plotcmd = ""
  sep = ""

  do for [it=1:words(TTOK)] {
    ttok = word(TTOK, it)
    tlbl_disp = fmt_time(word(TLBL_RAW, it))
    f = find_file(bh, ttok)

    if (strlen(f) == 0) {
      print sprintf("WARNING: Missing snapshot for overlay: bh=%s ttok=%s", bh, ttok)
    } else {
      # model (black line+points, varied dash)
      plotcmd = plotcmd . sep . \
        sprintf("'%s' every ::1 using (($2)*MODEL_SCALE):1 with linespoints ls %d title 'Model %s'", f, 20+it, tlbl_disp)

      # obs (red line+points)
      if (USE_TIDY_OBS) {
        system(make_obs_tmp(bh, ttok))
        nobs = obs_tmp_count(bh, ttok)
        if (nobs <= 0) {
          print sprintf("WARNING: No tidy OBS rows for overlay: bh=%s ttok=%s", bh, ttok)
        } else {
          plotcmd = plotcmd . ", " . \
            sprintf("'%s' using 1:2 with linespoints ls %d title 'Obs  %s'", obs_tmp_name(bh, ttok), 30+it, tlbl_disp)
        }
      } else {
        plotcmd = plotcmd . ", " . \
          sprintf("'%s' every ::1 using ($3):1 with linespoints ls %d title 'Obs  %s'", f, 30+it, tlbl_disp)
      }

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
set terminal pngcairo enhanced font "Arial,22" size 3000,650
unset key

do for [ib=1:words(BHs)] {
  bh = word(BHs, ib)

  # ---- compute a consistent x-range for this borehole (model only)
  xmin =  1e100
  xmax = -1e100
  have_any = 0

  do for [it=1:words(TTOK)] {
    ttok = word(TTOK, it)
    f = find_file(bh, ttok)
    if (strlen(f) > 0) {
      stats f every ::1 using (($2)*MODEL_SCALE) nooutput
      if (exists("STATS_valid") && STATS_valid > 0) {
        if (STATS_min < xmin) { xmin = STATS_min }
        if (STATS_max > xmax) { xmax = STATS_max }
        have_any = 1
      }
    }
  }

  if (have_any) {
    pad = 0.05*(xmax - xmin)
    if (pad <= 0) { pad = 1.0 }
    set xrange [xmin-pad : xmax+pad]
  } else {
    set xrange [*:*]
  }

  set format x "%.1f"
  set xtics autofreq

  set output sprintf("ERT_profile_%s_5panels.png", bh)
  set multiplot layout 1,5 rowsfirst

  do for [it=1:words(TTOK)] {
    ttok = word(TTOK, it)
    tlbl_disp = fmt_time(word(TLBL_RAW, it))
    f = find_file(bh, ttok)

    if (it == 1) {
      set ytics
      set ylabel "Depth (m)" font "Arial,26"
    } else {
      unset ytics
      set ylabel ""
    }

    set xlabel "Moisture" font "Arial,24"
    set title tlbl_disp font "Arial,24"

    if (strlen(f) == 0) {
      print sprintf("WARNING: Missing snapshot for panels: bh=%s ttok=%s", bh, ttok)
      eval(plot_missing(sprintf("MISSING\\n%s\\n%s", bh, ttok)))
    } else {
      unset label

      if (USE_TIDY_OBS) {
        system(make_obs_tmp(bh, ttok))
        nobs = obs_tmp_count(bh, ttok)

        if (nobs <= 0) {
          print sprintf("WARNING: No tidy OBS rows for panels: bh=%s ttok=%s", bh, ttok)
          plot \
            f every ::1 using (($2)*MODEL_SCALE):1 with linespoints ls 11 notitle
        } else {
          plot \
            f every ::1 using (($2)*MODEL_SCALE):1 with linespoints ls 11 notitle, \
            obs_tmp_name(bh, ttok) using 1:2        with linespoints ls 12 notitle
        }
      } else {
        plot \
          f every ::1 using (($2)*MODEL_SCALE):1 with linespoints ls 11 notitle, \
          f every ::1 using ($3):1               with linespoints ls 12 notitle
      }
    }
  }

  unset multiplot
  unset output
}

# ------------------------------------------------
# Cleanup temporary observation files
# ------------------------------------------------
system("rm -f __obs_*.dat")

unset output

