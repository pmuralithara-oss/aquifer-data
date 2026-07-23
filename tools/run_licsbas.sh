#!/bin/zsh
# Run the LiCSBAS chain on a prepared workdir containing GEOC/ (and optionally GACOS/).
# Usage: run_licsbas.sh <workdir> [n_para]
set -e
WORK="$1"
NP="${2:-8}"
REPO=/Users/pranavm/Desktop/aquifer-data/tools/LiCSBAS
export PATH="$HOME/micromamba/envs/licsbas/bin:$REPO/bin:$PATH"
export PYTHONPATH="$REPO/LiCSBAS_lib:$REPO/bin:$PYTHONPATH"
export MPLBACKEND=Agg
cd "$WORK"

echo "##### 02 ml_prep (nlook=1)"
LiCSBAS02_ml_prep.py -i GEOC -o GEOCml1 -n 1 --n_para "$NP"

IFGDIR=GEOCml1
if [ -d GACOS ] && [ "$(ls GACOS/*.sltd.geo.tif 2>/dev/null | wc -l | tr -d ' ')" -gt 0 ]; then
  echo "##### 03op GACOS"
  LiCSBAS03op_GACOS.py -i GEOCml1 -o GEOCml1GACOS -g GACOS --n_para "$NP" || {
    echo "GACOS step failed; continuing uncorrected"; }
  if [ -f GEOCml1GACOS/slc.mli.par ]; then IFGDIR=GEOCml1GACOS; else rm -rf GEOCml1GACOS; fi
else
  echo "##### no GACOS dir/files; skipping 03op"
fi
echo "##### using IFGDIR=$IFGDIR"

echo "##### 11 check_unw"
LiCSBAS11_check_unw.py -d "$IFGDIR"

echo "##### 12 loop_closure"
LiCSBAS12_loop_closure.py -d "$IFGDIR" --n_para "$NP"

echo "##### 13 sb_inv"
LiCSBAS13_sb_inv.py -d "$IFGDIR" --n_para "$NP" --mem_size 4096

TSDIR="TS_$IFGDIR"
echo "##### 14 vel_std"
LiCSBAS14_vel_std.py -t "$TSDIR" --mem_size 4096

echo "##### 15 mask_ts"
LiCSBAS15_mask_ts.py -t "$TSDIR"

echo "##### 16 filt_ts"
LiCSBAS16_filt_ts.py -t "$TSDIR" --n_para "$NP" --interpolate_nans

echo "##### DONE $WORK -> $TSDIR"
