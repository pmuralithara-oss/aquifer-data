#!/usr/bin/env python3
"""Render velocity / seasonal-amplitude / classification maps (PNG) from
decompose_ts.py GeoTIFF outputs, with town labels and scale bar."""
import argparse, csv, os
import numpy as np
from osgeo import gdal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm

gdal.UseExceptions()

def read(path):
    ds = gdal.Open(path)
    a = ds.GetRasterBand(1).ReadAsArray().astype(float)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    if nd is not None:
        a = np.where(a == nd, np.nan, a)
    gt = ds.GetGeoTransform()
    ext = [gt[0], gt[0] + gt[1] * ds.RasterXSize, gt[3] + gt[5] * ds.RasterYSize, gt[3]]
    return a, ext

def add_towns(ax, towns, ext):
    for nm, la, lo in towns:
        if ext[0] <= lo <= ext[1] and ext[2] <= la <= ext[3]:
            ax.plot(lo, la, "^", color="k", ms=5, mew=0.8, mfc="white")
            ax.annotate(nm, (lo, la), textcoords="offset points", xytext=(4, 4),
                        fontsize=7, weight="bold",
                        path_effects=None, color="k",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.65))

def add_scalebar(ax, ext, lat_mid):
    km = 5
    deg = km / (111.32 * np.cos(np.radians(lat_mid)))
    x0 = ext[0] + (ext[1] - ext[0]) * 0.06
    y0 = ext[2] + (ext[3] - ext[2]) * 0.06
    ax.plot([x0, x0 + deg], [y0, y0], "k-", lw=3, solid_capstyle="butt")
    ax.text(x0 + deg / 2, y0 + (ext[3] - ext[2]) * 0.015, f"{km} km",
            ha="center", fontsize=8, weight="bold")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="outputs path prefix e.g. outputs/AOI3_delhi/AOI3_asc027A")
    ap.add_argument("--towns", required=True, help="csv name,lat,lon")
    ap.add_argument("--title", required=True)
    ap.add_argument("--mli", help="optional background intensity tif (frame geo.mli clipped)")
    ap.add_argument("--vmax", type=float, default=None, help="symmetric vel color range (mm/yr)")
    a = ap.parse_args()

    towns = []
    for row in csv.DictReader(open(a.towns)):
        towns.append((row["name"], float(row["lat"]), float(row["lon"])))

    bg = None
    if a.mli and os.path.exists(a.mli):
        m, extm = read(a.mli)
        with np.errstate(invalid="ignore"):
            bg = np.log10(np.where(m > 0, m, np.nan))
        bglims = np.nanpercentile(bg, [2, 98])

    def render(tif, cmap, label, out, symmetric=False, vlim=None, categorical=False):
        arr, ext = read(tif)
        fig, ax = plt.subplots(figsize=(9.5, 7.2))
        if bg is not None:
            ax.imshow(bg, extent=extm, cmap="gray", vmin=bglims[0], vmax=bglims[1],
                      interpolation="nearest", aspect="auto")
        if categorical:
            cmap_c = ListedColormap(["#b2182b", "#2166ac", "#bdbdbd"])
            shown = np.where(arr == 0, np.nan, arr)
            im = ax.imshow(shown, extent=ext, cmap=cmap_c, vmin=1, vmax=3,
                           interpolation="nearest", aspect="auto", alpha=0.75)
            cb = fig.colorbar(im, ax=ax, shrink=0.75, ticks=[1.33, 2, 2.67])
            cb.ax.set_yticklabels(["inelastic\ncandidate", "elastic\ndominant", "indeterm."], fontsize=8)
        else:
            if symmetric:
                v = vlim or np.nanpercentile(np.abs(arr), 98)
                norm = TwoSlopeNorm(vcenter=0, vmin=-v, vmax=v)
                im = ax.imshow(arr, extent=ext, cmap=cmap, norm=norm,
                               interpolation="nearest", aspect="auto", alpha=0.8)
            else:
                lo, hi = np.nanpercentile(arr, [2, 98])
                im = ax.imshow(arr, extent=ext, cmap=cmap, vmin=lo, vmax=hi,
                               interpolation="nearest", aspect="auto", alpha=0.8)
            fig.colorbar(im, ax=ax, shrink=0.75, label=label)
        add_towns(ax, towns, ext)
        add_scalebar(ax, ext, (ext[2] + ext[3]) / 2)
        ax.set_title(a.title + " — " + label, fontsize=11)
        ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print("wrote", out)

    p = a.prefix
    if os.path.exists(f"{p}_los_rate_mmyr.tif"):
        render(f"{p}_los_rate_mmyr.tif", "RdBu_r", "LOS rate (mm/yr, + toward satellite)",
               f"{p}_los_rate_map.png", symmetric=True, vlim=a.vmax)
    if os.path.exists(f"{p}_quasi_vert_rate_mmyr.tif"):
        render(f"{p}_quasi_vert_rate_mmyr.tif", "RdBu_r", "quasi-vertical rate (mm/yr)",
               f"{p}_vert_rate_map.png", symmetric=True, vlim=a.vmax)
    if os.path.exists(f"{p}_seasonal_amp_mm.tif"):
        render(f"{p}_seasonal_amp_mm.tif", "viridis", "annual seasonal amplitude (mm)",
               f"{p}_seasonal_amp_map.png")
    if os.path.exists(f"{p}_class.tif"):
        render(f"{p}_class.tif", None, "elastic/inelastic classification",
               f"{p}_class_map.png", categorical=True)
    if os.path.exists(f"{p}_vert_rate_joint.tif"):
        render(f"{p}_vert_rate_joint.tif", "RdBu_r", "joint vertical rate (mm/yr, + up)",
               f"{p}_vert_rate_joint_map.png", symmetric=True, vlim=a.vmax)
    if os.path.exists(f"{p}_east_rate_joint.tif"):
        render(f"{p}_east_rate_joint.tif", "PuOr_r", "joint east rate (mm/yr, + east)",
               f"{p}_east_rate_joint_map.png", symmetric=True, vlim=a.vmax)

if __name__ == "__main__":
    main()
