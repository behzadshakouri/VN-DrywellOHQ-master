#!/usr/bin/gnuplot
# ============================================================
# Sensitivity overlay: MODEL ONLY across subfolders containing "sc"
#
# Outputs:
# (A) One plot per (borehole,time): overlays unique parameter sets
#     SENS_model_<BH>_<TLBL_FILE>.png
# (C) One plot per borehole: 5 subplots (1x5), each panel is a time,
#     each panel overlays unique parameter sets (model only)
#     SENS_model_<BH>_5panels.png
#
# Folder name parsing (tokens split by '_' or '-'):
#   sc <val>          -> global scale (only if no g/uw present)
#   g  <val>          -> gravel scale => sc_{g}=<val>
#   uw <val>          -> under-well  => sc_{uw}=<val>
#   nr <val>          -> default nr=12 if missing
#
# IMPORTANT: Your real format is like:
#   "53th sc g 3 uw 30 nr 16"
# which means g and uw are separate tokens (NOT "sc uw ...")
#
# NOTE: gnuplot cannot safely iterate directory names WITH SPACES
# because words() splits on whitespace. Use underscores/hyphens or symlinks.
#
# Snapshot CSV columns:
#   depth_m,theta_model,theta_obs,time
# We plot MODEL ONLY: x = theta_model * MODEL_SCALE, y = depth_m
# ============================================================

reset
set datafile separator ','
set grid
set tics out
set border lw 2

# Depth increases downward
set ylabel "Depth (m)" font "Arial,32"
set yrange [*:*] reverse
set xlabel "Soil Moisture" font "Arial,32"

# If model is 0..1 and you want percent, use 100; else 1.
MODEL_SCALE = 100

# Line widths
LW_DEFAULT = 3
LW_THICK   = 6   # for nr=16

# -------------------------------
# Inputs (edit if needed)
# -------------------------------
BHs  = "ERT-3 ERT-5"
TTOK = "45763p588194 45763p682639 45763p770833 45764p673611 45764p729167"

TLBL_RAW  = "2025-04-16_14:07 2025-04-16_16:23 2025-04-16_18:30 2025-04-17_16:10 2025-04-17_17:30"
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
fmt_time(s) = s[1:10] . " " . s[12:16]

# Find snapshot file in a given directory, try NEW then OLD
find_file_in_dir(dir, bh, ttok) = stripnl(system(sprintf( \
  "sh -c 'd=\"%s\"; " . \
  "f1=\"$d/%s%s%s%s%s\"; f2=\"$d/%s%s%s%s%s\"; " . \
  "if [ -f \"$f1\" ]; then echo \"$f1\"; " . \
  "elif [ -f \"$f2\" ]; then echo \"$f2\"; " . \
  "else echo \"\"; fi'", \
  dir, \
  PREFIX_NEW, ttok, MID_NEW, bh, SUFFIX, \
  PREFIX_OLD, ttok, MID_OLD, bh, SUFFIX)))

plot_missing(msg) = sprintf( \
  "unset label; set label 1 '%s' at graph 0.5, graph 0.5 center font 'Arial,28'; " . \
  "plot '-' using 1:2 with points pt 7 ps 0.01 notitle; 0 0; e", msg)

# List dirs containing "sc" (names must be underscore/hyphen separated; no spaces)
DIRS = stripnl(system("sh -c 'ls -1d *sc* 2>/dev/null | sort -V'"))

# Parse tokens from folder name (already normalized to spaces)
get_global_sc(tok) = stripnl(system(sprintf( \
  "sh -c 'echo \"%s\" | awk \"{for(i=1;i<=NF;i++) if(\\$i==\\\"sc\\\" && i<NF){if(\\$(i+1)!=\\\"g\\\" && \\$(i+1)!=\\\"uw\\\"){print \\$(i+1); exit}}}\"'", tok)))
get_g(tok) = stripnl(system(sprintf( \
  "sh -c 'echo \"%s\" | awk \"{for(i=1;i<=NF;i++) if(\\$i==\\\"g\\\" && i<NF){print \\$(i+1); exit}}\"'", tok)))
get_uw(tok) = stripnl(system(sprintf( \
  "sh -c 'echo \"%s\" | awk \"{for(i=1;i<=NF;i++) if(\\$i==\\\"uw\\\" && i<NF){print \\$(i+1); exit}}\"'", tok)))
get_nr(tok) = stripnl(system(sprintf( \
  "sh -c 'echo \"%s\" | awk \"{for(i=1;i<=NF;i++) if(\\$i==\\\"nr\\\" && i<NF){print \\$(i+1); exit}}\"'", tok)))

# Build label with subscripts (enhanced text) from parsed values
# - if no g and no uw -> global sc
# - else use sc_{g}, sc_{uw}
make_label(sc, scg, scuw, nr) = stripnl(system(sprintf( \
  "sh -c 'sc=\"%s\"; scg=\"%s\"; scuw=\"%s\"; nr=\"%s\"; " . \
  "if [ -z \"$nr\" ]; then nr=12; fi; " . \
  "lab=\"\"; " . \
  "if [ -z \"$scg\" ] && [ -z \"$scuw\" ]; then " . \
  "  [ -z \"$sc\" ] && sc=\"?\"; " . \
  "  lab=\"sc=$sc (global)\"; " . \
  "else " . \
  "  if [ -n \"$scg\" ]; then lab=\"${lab}sc_{g}=$scg \"; fi; " . \
  "  if [ -n \"$scuw\" ]; then lab=\"${lab}sc_{uw}=$scuw \"; fi; " . \
  "fi; " . \
  "if [ \"$nr\" != \"12\" ]; then lab=\"${lab}nr=$nr\"; fi; " . \
  "echo \"$lab\" | sed \"s/[[:space:]]*$//\"' ", \
  sc, scg, scuw, nr)))

# ============================================================
# (A) One plot per (borehole,time): overlay unique scale combos
# ============================================================
set terminal pngcairo enhanced font "Arial,26" size 1600,900
set key outside right top
set key font "Arial,20"
set key samplen 2 spacing 1.1
set key opaque box

if (strlen(DIRS) == 0) {
  print "ERROR: No folders found matching *sc* in current directory."
  set output "SENS_NO_FOLDERS.png"
  eval(plot_missing("No *sc* folders found"))
  unset output
} else {

  array USED_LABELS_A[400]

  do for [ib=1:words(BHs)] {
    bh = word(BHs, ib)

    do for [it=1:words(TTOK)] {
      ttok = word(TTOK, it)
      tlbl_disp = fmt_time(word(TLBL_RAW, it))
      tlbl_file = word(TLBL_FILE, it)

      set output sprintf("SENS_model_%s_%s.png", bh, tlbl_file)
      set title sprintf("%s | %s | MODEL sensitivity", bh, tlbl_disp) font "Arial,30"

      plotcmd = ""
      sep = ""
      n_used = 0

      do for [id=1:words(DIRS)] {
        d = word(DIRS, id)
        f = find_file_in_dir(d, bh, ttok)
        if (strlen(f) == 0) { continue }

        tok = stripnl(system(sprintf("sh -c 'echo \"%s\" | tr \"_-\" \"  \"'", d)))

        sc   = get_global_sc(tok)
        scg  = get_g(tok)
        scuw = get_uw(tok)
        nr   = get_nr(tok)
        if (strlen(nr) == 0) { nr = "12" }

        lab = make_label(sc, scg, scuw, nr)

        # deduplicate by label
        duplicate = 0
        do for [k=1:n_used] {
          if (USED_LABELS_A[k] eq lab) { duplicate = 1 }
        }
        if (duplicate) { continue }

        n_used = n_used + 1
        USED_LABELS_A[n_used] = lab

        dtv = 1 + (n_used-1) % 10
        lwv = LW_DEFAULT
        if (nr == "16") { lwv = LW_THICK }

        plotcmd = plotcmd . sep . \
          sprintf("'%s' every ::1 using (($2)*MODEL_SCALE):1 with lines lw %d dt %d title '%s'", \
                  f, lwv, dtv, lab)
        sep = ", "
      }

      if (strlen(plotcmd) > 0) {
        unset label
        eval("plot " . plotcmd)
      } else {
        eval(plot_missing(sprintf("NO SNAPSHOTS FOUND\\n%s\\n%s", bh, ttok)))
      }

      unset output
    }
  }
}

# ============================================================
# (C) Per borehole: 5 subplots (1x5) — NO LEGEND
#   - Keeps all panels equal size; avoids overlap/clutter.
#   - Still uses thicker lines for nr=16.
# ============================================================
set terminal pngcairo enhanced font "Arial,22" size 3400,750
unset key

if (strlen(DIRS) > 0) {

  array USED_LABELS_C[400]

  do for [ib=1:words(BHs)] {
    bh = word(BHs, ib)

    set output sprintf("SENS_model_%s_5panels.png", bh)
    set multiplot layout 1,5 rowsfirst

    do for [it=1:words(TTOK)] {

      ttok = word(TTOK, it)
      tlbl_disp = fmt_time(word(TLBL_RAW, it))

      # y-axis only on first panel
      if (it == 1) {
        set ylabel "Depth (m)" font "Arial,26"
        set ytics
      } else {
        unset ytics
        set ylabel ""
      }

      set xlabel "Moisture" font "Arial,24"
      set title tlbl_disp font "Arial,24"

      plotcmd = ""
      sep = ""
      n_used = 0

      do for [id=1:words(DIRS)] {
        d = word(DIRS, id)
        f = find_file_in_dir(d, bh, ttok)
        if (strlen(f) == 0) { continue }

        tok = stripnl(system(sprintf("sh -c 'echo \"%s\" | tr \"_-\" \"  \"'", d)))

        sc   = get_global_sc(tok)
        scg  = get_g(tok)
        scuw = get_uw(tok)
        nr   = get_nr(tok)
        if (strlen(nr) == 0) { nr = "12" }

        lab = make_label(sc, scg, scuw, nr)

        # deduplicate by label (per panel, just to keep dt index stable-ish)
        duplicate = 0
        do for [k=1:n_used] {
          if (USED_LABELS_C[k] eq lab) { duplicate = 1 }
        }
        if (duplicate) { continue }

        n_used = n_used + 1
        USED_LABELS_C[n_used] = lab

        dtv = 1 + (n_used-1) % 10
        lwv = LW_DEFAULT
        if (nr == "16") { lwv = LW_THICK }

        plotcmd = plotcmd . sep . \
          sprintf("'%s' every ::1 using (($2)*MODEL_SCALE):1 with lines lw %d dt %d notitle", \
                  f, lwv, dtv)

        sep = ", "
      }

      if (strlen(plotcmd) > 0) {
        unset label
        eval("plot " . plotcmd)
      } else {
        eval(plot_missing(sprintf("NO SNAPSHOTS\\n%s\\n%s", bh, ttok)))
      }
    }

    unset multiplot
    unset output
  }
}

unset output
