#!/usr/bin/env python3
"""Align all pair tifs in a GEOC dir onto the majority grid (fixes LiCSAR sub-pixel
grid shifts between processing eras). Mismatched files are warped (nearest) in place."""
import argparse, glob, os
from collections import Counter
from osgeo import gdal

gdal.UseExceptions()

def grid_of(path):
    ds = gdal.Open(path)
    gt = ds.GetGeoTransform()
    return (ds.RasterXSize, ds.RasterYSize, round(gt[0], 6), round(gt[3], 6),
            round(gt[1], 6), round(gt[5], 6))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geoc")
    a = ap.parse_args()
    tifs = sorted(glob.glob(os.path.join(a.geoc, "20*_20*", "*.geo.*.tif")))
    grids = Counter(grid_of(t) for t in tifs)
    ref = grids.most_common(1)[0][0]
    W, H, x0, y0, dx, dy = ref
    bounds = (x0, y0 + H * dy, x0 + W * dx, y0)  # (minx, miny, maxx, maxy)
    print(f"majority grid {W}x{H} origin ({x0},{y0}); {len(tifs)} files, "
          f"{sum(c for g, c in grids.items() if g != ref)} to align")
    fixed = 0
    for t in tifs:
        if grid_of(t) == ref:
            continue
        tmp = t + ".aligned.tif"
        gdal.Warp(tmp, t, outputBounds=bounds, xRes=dx, yRes=abs(dy),
                  resampleAlg="near", creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])
        os.replace(tmp, t)
        fixed += 1
    print(f"aligned {fixed} files")
    # frame-level tifs too
    for t in sorted(glob.glob(os.path.join(a.geoc, "*.geo.*.tif"))):
        if grid_of(t) != ref:
            tmp = t + ".aligned.tif"
            gdal.Warp(tmp, t, outputBounds=bounds, xRes=dx, yRes=abs(dy),
                      resampleAlg="near", creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])
            os.replace(tmp, t)
            print(f"aligned frame file {os.path.basename(t)}")

if __name__ == "__main__":
    main()
