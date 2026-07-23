#!/usr/bin/env python3
"""Crop a clipped GEOC/GACOS tree to a smaller AOI box (local, no network)."""
import argparse, glob, os, shutil
import numpy as np
from osgeo import gdal

gdal.UseExceptions()

def crop_tif(src, dst, box):
    lat_n, lat_s, lon_w, lon_e = box
    ds = gdal.Open(src)
    gdal.Translate(dst, ds, projWin=[lon_w, lat_n, lon_e, lat_s],
                   creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])
    ds = None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, help="dir containing GEOC/ (and GACOS/)")
    ap.add_argument("--dst", required=True)
    ap.add_argument("--box", required=True, help="latN,latS,lonW,lonE")
    a = ap.parse_args()
    box = [float(x) for x in a.box.split(",")]
    sg, dg = os.path.join(a.src, "GEOC"), os.path.join(a.dst, "GEOC")
    os.makedirs(dg, exist_ok=True)
    # frame-level files
    for f in glob.glob(os.path.join(sg, "*")):
        b = os.path.basename(f)
        if os.path.isdir(f):
            continue
        if f.endswith(".tif"):
            crop_tif(f, os.path.join(dg, b), box)
        else:
            shutil.copy2(f, os.path.join(dg, b))
    n = 0
    for pdir in sorted(glob.glob(os.path.join(sg, "20*_20*"))):
        pr = os.path.basename(pdir)
        od = os.path.join(dg, pr)
        os.makedirs(od, exist_ok=True)
        ok = True
        for f in glob.glob(os.path.join(pdir, "*.tif")):
            try:
                crop_tif(f, os.path.join(od, os.path.basename(f)), box)
            except Exception as e:
                print(f"WARN {pr}: {e}"); ok = False
        if not ok:
            shutil.rmtree(od)
        else:
            n += 1
        if n % 200 == 0:
            print(f"  {n} pairs cropped", flush=True)
    print(f"cropped {n} pairs -> {dg}")
    sgac = os.path.join(a.src, "GACOS")
    if os.path.isdir(sgac):
        dgac = os.path.join(a.dst, "GACOS")
        os.makedirs(dgac, exist_ok=True)
        m = 0
        for f in sorted(glob.glob(os.path.join(sgac, "*.sltd.geo.tif"))):
            try:
                crop_tif(f, os.path.join(dgac, os.path.basename(f)), box); m += 1
            except Exception as e:
                print(f"WARN gacos {os.path.basename(f)}: {e}")
        print(f"cropped {m} GACOS epochs")

if __name__ == "__main__":
    main()
