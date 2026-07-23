#!/usr/bin/env python3
"""Joint asc+desc LOS -> (vertical, east) velocity decomposition on the common
time window and common 0.001-deg grid. Both TS dirs must come from the same clip box.

Per pixel: [rate_asc]   [U_a  E_a] [v_up ]
           [rate_desc] = [U_d  E_d] [v_east]   (north term neglected: S1 insensitivity)
"""
import argparse, os
import numpy as np
import h5py
from osgeo import gdal, osr
import matplotlib
matplotlib.use("Agg")

gdal.UseExceptions()

def load(tsdir, cumfile="cum.h5"):
    with h5py.File(os.path.join(tsdir, cumfile), "r") as h:
        d = {
            "cum": h["cum"][()].astype(np.float32),
            "dates": [str(int(x)) for x in h["imdates"][()]],
            "clat": float(h["corner_lat"][()]), "clon": float(h["corner_lon"][()]),
            "plat": float(h["post_lat"][()]), "plon": float(h["post_lon"][()]),
            "E": h["E.geo"][()], "N": h["N.geo"][()], "U": h["U.geo"][()],
        }
    m = os.path.join(tsdir, "results", "mask")
    sh = d["cum"].shape[1:]
    d["mask"] = np.fromfile(m, np.float32).reshape(sh) if os.path.exists(m) else np.ones(sh, np.float32)
    return d

def dec_year(dates):
    from datetime import datetime as dt
    out = []
    for s in dates:
        t = dt.strptime(s, "%Y%m%d")
        y0, y1 = dt(t.year, 1, 1), dt(t.year + 1, 1, 1)
        out.append(t.year + (t - y0).total_seconds() / (y1 - y0).total_seconds())
    return np.array(out)

def fit_rate(t, Y):
    G = np.stack([np.ones_like(t), t, np.sin(2*np.pi*t), np.cos(2*np.pi*t),
                  np.sin(4*np.pi*t), np.cos(4*np.pi*t)], 1)
    n_px = Y.shape[1]
    rate = np.full(n_px, np.nan, np.float32)
    CH = 20000
    for s in range(0, n_px, CH):
        e = min(n_px, s + CH)
        Yc = Y[:, s:e]; M = np.isfinite(Yc)
        Y0 = np.where(M, Yc, 0.0)
        GtG = np.einsum("ep,eq,en->pqn", G, G, M.astype(np.float64), optimize=True)
        Gty = np.einsum("ep,en->pn", G, Y0, optimize=True)
        nob = M.sum(0)
        for i in range(e - s):
            if nob[i] < 16: continue
            try: rate[s + i] = np.linalg.solve(GtG[:, :, i], Gty[:, i])[1]
            except np.linalg.LinAlgError: pass
    return rate

def geotiff(path, arr, gt):
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32, options=["COMPRESS=DEFLATE"])
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference(); srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(arr); ds.FlushCache()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asc", required=True); ap.add_argument("--desc", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--refbox", help="latN,latS,lonW,lonE common stable-zone re-reference")
    ap.add_argument("--deramp", action="store_true", help="remove planar ramp from each rate field before solving")
    a = ap.parse_args()
    A, D = load(a.asc), load(a.desc)
    if a.refbox:
        rn, rs, rw, re_ = [float(x) for x in a.refbox.split(",")]
        for X in (A, D):
            n_im = X["cum"].shape[0]
            i0 = max(0, int((X["clat"] - rn) / -X["plat"])); i1 = int((X["clat"] - rs) / -X["plat"])
            j0 = max(0, int((rw - X["clon"]) / X["plon"])); j1 = int((re_ - X["clon"]) / X["plon"])
            refs = np.nanmedian(X["cum"][:, i0:i1, j0:j1].reshape(n_im, -1), axis=1)
            X["cum"] = X["cum"] - refs[:, None, None]
        print(f"re-referenced both geometries to {a.refbox}")
    # align D onto A's lattice by integer pixel shift (same 0.001-deg spacing)
    dy = int(round((A["clat"] - D["clat"]) / -A["plat"]))   # +ve: D starts dy rows above A
    dx = int(round((A["clon"] - D["clon"]) / A["plon"]))
    HA, WA = A["cum"].shape[1:]; HD, WD = D["cum"].shape[1:]
    # overlap in A-coordinates
    ya0 = max(0, -dy); xa0 = max(0, -dx)
    ya1 = min(HA, HD - dy); xa1 = min(WA, WD - dx)
    if ya1 <= ya0 or xa1 <= xa0:
        raise SystemExit("no overlap between asc and desc grids")
    def crop_A(x): return x[..., ya0:ya1, xa0:xa1]
    def crop_D(x): return x[..., ya0 + dy:ya1 + dy, xa0 + dx:xa1 + dx]
    for k in ("cum", "E", "N", "U", "mask"):
        A[k] = crop_A(A[k]); D[k] = crop_D(D[k])
    A["clat"] += ya0 * A["plat"]; A["clon"] += xa0 * A["plon"]
    print(f"grid align: shift dy={dy} dx={dx}, overlap {ya1-ya0}x{xa1-xa0}")
    sh = A["cum"].shape[1:]
    tA, tD = dec_year(A["dates"]), dec_year(D["dates"])
    lo = max(tA.min(), tD.min()); hi = min(tA.max(), tD.max())
    sA = (tA >= lo) & (tA <= hi); sD = (tD >= lo) & (tD <= hi)
    print(f"common window {lo:.2f}-{hi:.2f}: asc {sA.sum()} epochs, desc {sD.sum()} epochs")
    YA = np.where(A["mask"].reshape(-1) == 1, A["cum"][sA].reshape(sA.sum(), -1), np.nan)
    YD = np.where(D["mask"].reshape(-1) == 1, D["cum"][sD].reshape(sD.sum(), -1), np.nan)
    rA = fit_rate(tA[sA] - lo, YA)
    rD = fit_rate(tD[sD] - lo, YD)
    if a.deramp:
        H, W = sh
        yy, xx = np.mgrid[0:H, 0:W]
        for r in (rA, rD):
            rr = r.reshape(H, W)
            m = np.isfinite(rr)
            Gd = np.stack([np.ones(m.sum()), xx[m], yy[m]], 1)
            c = np.linalg.lstsq(Gd, rr[m], rcond=None)[0]
            rr -= (c[0] + c[1] * xx + c[2] * yy)
        print("deramped both rate fields (planar)")
    Ua, Ea = A["U"].reshape(-1), A["E"].reshape(-1)
    Ud, Ed = D["U"].reshape(-1), D["E"].reshape(-1)
    det = Ua * Ed - Ea * Ud
    ok = np.isfinite(rA) & np.isfinite(rD) & (np.abs(det) > 0.05)
    vU = np.full(rA.shape, np.nan, np.float32); vE = np.full(rA.shape, np.nan, np.float32)
    vU[ok] = (Ed[ok] * rA[ok] - Ea[ok] * rD[ok]) / det[ok]
    vE[ok] = (-Ud[ok] * rA[ok] + Ua[ok] * rD[ok]) / det[ok]
    gt = (A["clon"], A["plon"], 0.0, A["clat"], 0.0, A["plat"])
    os.makedirs(a.out, exist_ok=True)
    geotiff(os.path.join(a.out, f"{a.name}_vert_rate_joint.tif"), vU.reshape(sh), gt)
    geotiff(os.path.join(a.out, f"{a.name}_east_rate_joint.tif"), vE.reshape(sh), gt)
    print(f"joint px: {ok.sum()}; vert p2/p50/p98 = "
          f"{np.nanpercentile(vU[ok],2):.1f}/{np.nanpercentile(vU[ok],50):.1f}/{np.nanpercentile(vU[ok],98):.1f} mm/yr")

if __name__ == "__main__":
    main()
