#!/usr/bin/env python3
"""Integrity-check a clipped GEOC stack: every pair's unw+cc must open, match the
reference grid, and contain plausible data. Corrupt/suspect files are deleted
(so a re-run of fetch_ceda_clip.py refetches them) and listed."""
import argparse, glob, os, sys
import numpy as np
from osgeo import gdal

gdal.UseExceptions()

def check(path, refshape):
    try:
        ds = gdal.Open(path)
        if (ds.RasterYSize, ds.RasterXSize) != refshape:
            return f"shape {(ds.RasterYSize, ds.RasterXSize)} != {refshape}"
        a = ds.GetRasterBand(1).ReadAsArray()
        if a is None:
            return "read failed"
        fin = np.isfinite(a)
        if fin.sum() == 0:
            return "all non-finite"
        nz = (a[fin] != 0).mean() if fin.sum() else 0
        if nz < 0.01:
            return f"only {nz*100:.2f}% nonzero"
        return None
    except Exception as e:
        return f"exception: {e}"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geoc")
    ap.add_argument("--delete", action="store_true")
    a = ap.parse_args()
    # reference grid from frame E tif
    etif = glob.glob(os.path.join(a.geoc, "*.geo.E.tif"))
    if not etif:
        sys.exit("no frame E tif found")
    ds = gdal.Open(etif[0]); refshape = (ds.RasterYSize, ds.RasterXSize)
    bad = []
    pairs = sorted(glob.glob(os.path.join(a.geoc, "20*_20*")))
    for p in pairs:
        pr = os.path.basename(p)
        for kind in ("geo.unw.tif", "geo.cc.tif"):
            f = os.path.join(p, f"{pr}.{kind}")
            if not os.path.exists(f):
                bad.append((f, "missing")); continue
            err = check(f, refshape)
            if err:
                bad.append((f, err))
    print(f"{len(pairs)} pairs checked against grid {refshape}; {len(bad)} bad files")
    for f, e in bad[:40]:
        print("  BAD", os.path.basename(f), "->", e)
    if a.delete:
        for f, e in bad:
            if os.path.exists(f):
                os.remove(f)
        print(f"deleted {sum(1 for f,_ in bad if not os.path.exists(f))} files (refetch with downloader)")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
