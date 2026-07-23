#!/usr/bin/env python3
"""
Phase-2 decomposition of LiCSBAS time series: separate secular (candidate-inelastic)
from seasonal (elastic) deformation per pixel, with bootstrap uncertainties and a
classification per the mission spec.

Model per pixel (t in decimal years since first epoch):
    d(t) = a + b*t + c1*sin(2πt) + c2*cos(2πt) [+ s1*sin(4πt) + s2*cos(4πt)]
Semiannual terms kept only where they improve BIC.

Classification:
  inelastic_candidate: b ≤ -10 mm/yr AND |b| > 2σ_b AND rate in each half of the series
                       ≤ -5 mm/yr (sustained) AND ≥3 post-monsoon (Oct-Dec) yearly maxima
                       with ≥75% of successive steps downward.
  elastic_dominant:    annual amplitude ≥ 2×|b×1yr| AND |b| ≤ 2 mm/yr.
  indeterminate:       everything else (incl. uplift, noisy).
  (classified only where the LiCSBAS mask is valid and enough epochs exist)

Outputs (per run): GeoTIFFs (LOS rate, σ, annual amp, phase, class, R²),
decomposition_results.csv, hotspot time-series plots, summary json.
"""
import argparse, json, os, sys
import numpy as np
import h5py
from osgeo import gdal, osr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime as dt

gdal.UseExceptions()

def dec_year(dates):
    out = []
    for d in dates:
        d = dt.strptime(str(d), "%Y%m%d")
        y0 = dt(d.year, 1, 1)
        y1 = dt(d.year + 1, 1, 1)
        out.append(d.year + (d - y0).total_seconds() / (y1 - y0).total_seconds())
    return np.array(out)

def build_G(t, semi):
    cols = [np.ones_like(t), t, np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)]
    if semi:
        cols += [np.sin(4 * np.pi * t), np.cos(4 * np.pi * t)]
    return np.stack(cols, axis=1)

def fit_all(t, Y, semi):
    """Y: (n_ep, n_px) with NaNs. Returns params (n_par, n_px), resid, dof, rms."""
    G = build_G(t, semi)
    n_par = G.shape[1]
    n_ep, n_px = Y.shape
    params = np.full((n_par, n_px), np.nan, np.float64)
    rms = np.full(n_px, np.nan)
    r2 = np.full(n_px, np.nan)
    finite = np.isfinite(Y)
    nobs = finite.sum(0)
    # group pixels by identical availability pattern would be ideal; pragmatic: loop over
    # unique patterns is heavy — instead fill NaNs with model-free means using masked lstsq
    # via per-pixel normal equations vectorised in chunks.
    CH = 20000
    for s in range(0, n_px, CH):
        e = min(n_px, s + CH)
        Yc = Y[:, s:e]
        M = np.isfinite(Yc)
        Y0 = np.where(M, Yc, 0.0)
        # normal equations with per-pixel masking: A = G^T W G, b = G^T W y
        # G: (n_ep,n_par); M: (n_ep,npx)
        GtG = np.einsum("ep,eq,en->pqn", G, G, M.astype(np.float64), optimize=True)
        Gty = np.einsum("ep,en->pn", G, Y0, optimize=True)
        for i in range(e - s):
            if nobs[s + i] < n_par + 10:
                continue
            try:
                params[:, s + i] = np.linalg.solve(GtG[:, :, i], Gty[:, i])
            except np.linalg.LinAlgError:
                continue
        fit = G @ params[:, s:e]
        res = np.where(M, Yc - fit, np.nan)
        rms[s:e] = np.sqrt(np.nanmean(res ** 2, 0))
        var = np.nanvar(Yc, 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            r2[s:e] = 1.0 - np.nanmean(res ** 2, 0) / np.where(var > 0, var, np.nan)
    return params, rms, r2, nobs

def bic(rms, nobs, n_par):
    with np.errstate(invalid="ignore", divide="ignore"):
        return nobs * np.log(rms ** 2) + n_par * np.log(nobs)

def block_bootstrap(t, Y, params, semi, n_boot=150, block_days=90, seed=42):
    """σ of trend + annual amplitude via temporal block bootstrap of residuals."""
    rng = np.random.default_rng(seed)
    G = build_G(t, semi)
    fit = G @ params
    res = Y - fit
    n_ep, n_px = Y.shape
    # temporal blocks
    tdays = (t - t[0]) * 365.25
    edges = np.arange(0, tdays[-1] + block_days, block_days)
    blocks = [np.where((tdays >= a) & (tdays < b))[0] for a, b in zip(edges[:-1], edges[1:])]
    blocks = [b for b in blocks if len(b) > 0]
    Ginv = np.linalg.pinv(G)  # (n_par, n_ep)
    bs_b = np.empty((n_boot, n_px), np.float32)
    bs_amp = np.empty((n_boot, n_px), np.float32)
    for k in range(n_boot):
        order = rng.integers(0, len(blocks), len(blocks))
        idx = np.concatenate([blocks[i] for i in order])[:n_ep]
        if len(idx) < n_ep:
            idx = np.concatenate([idx, rng.integers(0, n_ep, n_ep - len(idx))])
        res_k = res[idx, :]
        # re-fit on original times with resampled residuals added to model
        Yk = fit + np.where(np.isfinite(res_k), res_k, 0.0)
        mk = Ginv @ np.where(np.isfinite(Yk), Yk, 0.0)
        bs_b[k] = mk[1]
        bs_amp[k] = np.sqrt(mk[2] ** 2 + mk[3] ** 2)
    return np.nanstd(bs_b, 0), np.nanstd(bs_amp, 0)

def post_monsoon_test(t, dates_dt, Y, min_years=3):
    """Fraction of successive Oct-Dec yearly maxima that step downward.
    Returns (frac_down, n_years). Y (n_ep, n_px)."""
    years = np.array([d.year for d in dates_dt])
    months = np.array([d.month for d in dates_dt])
    sel_pm = (months >= 10) & (months <= 12)
    yr_list = sorted(set(years[sel_pm]))
    maxima = []
    for y in yr_list:
        sel = sel_pm & (years == y)
        if sel.sum() == 0:
            continue
        maxima.append(np.nanmax(Y[sel, :], axis=0))
    if len(maxima) < min_years:
        return None, len(maxima)
    Mx = np.stack(maxima)  # (n_yr, n_px)
    steps = np.diff(Mx, axis=0)
    with np.errstate(invalid="ignore"):
        frac_down = np.nanmean((steps < 0).astype(float), axis=0)
    return frac_down, Mx.shape[0]

def geotiff(path, arr, gt, nodata=np.nan, dtype=gdal.GDT_Float32):
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, arr.shape[1], arr.shape[0], 1, dtype,
                    options=["COMPRESS=DEFLATE"])
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference(); srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1)
    if nodata is not None and not (isinstance(nodata, float) and np.isnan(nodata)):
        b.SetNoDataValue(float(nodata))
    b.WriteArray(np.asarray(arr))
    ds.FlushCache()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ts", required=True, help="LiCSBAS TS dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", required=True, help="label e.g. AOI3_delhi_asc")
    ap.add_argument("--cumfile", default="cum.h5", help="cum.h5 (unfiltered) or cum_filt.h5")
    ap.add_argument("--hotspots", help="csv: name,lat,lon[,published_rate_mm]")
    ap.add_argument("--n_boot", type=int, default=150)
    ap.add_argument("--utif", help="frame geo.U.tif (clipped) for LOS->quasi-vertical")
    ap.add_argument("--refbox", help="latN,latS,lonW,lonE: re-reference cum to per-epoch "
                    "median of this box (use a known-stable zone)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    with h5py.File(os.path.join(a.ts, a.cumfile), "r") as h:
        cum = h["cum"][()].astype(np.float32)  # mm, (n_im, len, wid)
        imdates = [str(int(x)) for x in h["imdate"][()]] if "imdate" in h else \
                  [str(int(x)) for x in h["imdates"][()]]
        vel_ref = h["vel"][()] if "vel" in h else None
        corner_lat = float(h["corner_lat"][()]); corner_lon = float(h["corner_lon"][()])
        post_lat = float(h["post_lat"][()]); post_lon = float(h["post_lon"][()])
        Ugeo = h["U.geo"][()] if "U.geo" in h else None
    n_im, length, width = cum.shape
    gt = (corner_lon, post_lon, 0.0, corner_lat, 0.0, post_lat)

    if a.refbox:
        rn, rs, rw, re_ = [float(x) for x in a.refbox.split(",")]
        i0 = max(0, int((corner_lat - rn) / -post_lat)); i1 = min(length, int((corner_lat - rs) / -post_lat))
        j0 = max(0, int((rw - corner_lon) / post_lon)); j1 = min(width, int((re_ - corner_lon) / post_lon))
        refseries = np.nanmedian(cum[:, i0:i1, j0:j1].reshape(n_im, -1), axis=1)
        nref = np.isfinite(cum[0, i0:i1, j0:j1]).sum()
        print(f"re-referencing to box {a.refbox} (px window {j0}:{j1}/{i0}:{i1}, {nref} px): "
              f"ref-series trend removed")
        cum = cum - refseries[:, None, None]

    maskf = os.path.join(a.ts, "results", "mask")
    mask = np.fromfile(maskf, np.float32).reshape(length, width) if os.path.exists(maskf) else np.ones((length, width), np.float32)
    cohf = os.path.join(a.ts, "results", "coh_avg")
    coh = np.fromfile(cohf, np.float32).reshape(length, width) if os.path.exists(cohf) else None

    dates_dt = [dt.strptime(d, "%Y%m%d") for d in imdates]
    t = dec_year(imdates); t = t - t[0]
    span = t[-1]
    Y = cum.reshape(n_im, -1)
    ok_mask = (mask.reshape(-1) == 1) & (np.isfinite(Y).sum(0) >= max(30, n_im * 0.4))
    print(f"{a.name}: {n_im} epochs over {span:.2f} yr; {ok_mask.sum()} / {length*width} pixels usable")
    Yv = np.where(ok_mask[None, :], Y, np.nan)

    # fit without and with semiannual; keep semiannual where BIC improves
    p4, rms4, r2_4, nobs = fit_all(t, Yv, semi=False)
    p6, rms6, r2_6, _ = fit_all(t, Yv, semi=True)
    use6 = bic(rms6, nobs, 6) < bic(rms4, nobs, 4)
    print(f"semiannual retained on {np.nansum(use6 & ok_mask)/max(1,ok_mask.sum())*100:.1f}% of usable pixels")

    b_rate = np.where(use6, p6[1], p4[1])
    amp_ann = np.where(use6, np.hypot(p6[2], p6[3]), np.hypot(p4[2], p4[3]))
    # phase of annual peak (decimal month of maximum of c1 sin + c2 cos)
    c1 = np.where(use6, p6[2], p4[2]); c2 = np.where(use6, p6[3], p4[3])
    phase_peak_yrfrac = (np.arctan2(c2, c1) / (2 * np.pi)) % 1.0  # sin(2πt+φ) peak at t=(0.25-φ/2π)
    # peak time: maximize c1 sin + c2 cos = A sin(2πt + φ0), φ0=atan2(c2,c1); peak t = (0.25 - φ0/(2π)) mod 1
    phase_peak_yrfrac = (0.25 - np.arctan2(c2, c1) / (2 * np.pi)) % 1.0
    rms = np.where(use6, rms6, rms4)
    r2 = np.where(use6, r2_6, r2_4)

    # bootstrap σ (use the 6-par model for σ everywhere for simplicity/conservatism)
    sig_b, sig_amp = block_bootstrap(t, np.where(np.isfinite(Yv), Yv, np.nan), p6, True,
                                     n_boot=a.n_boot)

    # sustained test: rate in each half
    half = t < (span / 2)
    pA, _, _, _ = fit_all(t[half], Yv[half], semi=True)
    pB, _, _, _ = fit_all(t[~half], Yv[~half], semi=True)
    rateA, rateB = pA[1], pB[1]

    fd, nyr = post_monsoon_test(t, dates_dt, Yv)
    frac_down = fd if fd is not None else np.full(length * width, np.nan)

    CLS = np.zeros(length * width, np.int16)  # 0 nodata,1 inelastic,2 elastic,3 indeterminate
    valid = ok_mask & np.isfinite(b_rate)
    inel = valid & (b_rate <= -10) & (np.abs(b_rate) > 2 * sig_b) & \
           (rateA <= -5) & (rateB <= -5) & (frac_down >= 0.75) & (nyr >= 3)
    elas = valid & ~inel & (amp_ann >= 2 * np.abs(b_rate)) & (np.abs(b_rate) <= 2)
    CLS[valid] = 3
    CLS[elas] = 2
    CLS[inel] = 1
    print(f"classified: inelastic {inel.sum()}, elastic {elas.sum()}, indeterminate {int(valid.sum()-inel.sum()-elas.sum())}")

    sh = (length, width)
    out = lambda n: os.path.join(a.out, f"{a.name}_{n}")
    geotiff(out("los_rate_mmyr.tif"), b_rate.reshape(sh), gt)
    geotiff(out("los_rate_sigma.tif"), sig_b.reshape(sh), gt)
    geotiff(out("seasonal_amp_mm.tif"), amp_ann.reshape(sh), gt)
    geotiff(out("seasonal_peak_yrfrac.tif"), phase_peak_yrfrac.reshape(sh), gt)
    geotiff(out("fit_r2.tif"), r2.reshape(sh), gt)
    geotiff(out("class.tif"), CLS.reshape(sh).astype(np.int16), gt, nodata=0, dtype=gdal.GDT_Int16)

    # quasi-vertical
    u = None
    if Ugeo is not None and Ugeo.shape == sh:
        u = Ugeo
    elif a.utif and os.path.exists(a.utif):
        cand = gdal.Open(a.utif).ReadAsArray()
        if cand.shape == sh:
            u = cand
    if u is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            geotiff(out("quasi_vert_rate_mmyr.tif"),
                    (b_rate.reshape(sh) / np.where(np.abs(u) > 0.3, u, np.nan)), gt)

    # CSV (valid pixels only)
    ii, jj = np.divmod(np.where(valid)[0], width)
    lats = corner_lat + (ii + 0.5) * post_lat
    lons = corner_lon + (jj + 0.5) * post_lon
    cls_names = np.array(["nodata", "inelastic_candidate", "elastic_dominant", "indeterminate"])
    import csv
    with open(out("decomposition_results.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pid", "lat", "lon", "los_rate_mmyr", "rate_sigma", "seasonal_amp_mm",
                    "amp_sigma", "peak_yrfrac", "r2", "rms_mm", "coh_avg",
                    "rate_1sthalf", "rate_2ndhalf", "frac_postmonsoon_down", "class"])
        k = np.where(valid)[0]
        for n, (pid, la, lo) in enumerate(zip(k, lats, lons)):
            w.writerow([int(pid), f"{la:.5f}", f"{lo:.5f}", f"{b_rate[pid]:.2f}",
                        f"{sig_b[pid]:.2f}", f"{amp_ann[pid]:.2f}", f"{sig_amp[pid]:.2f}",
                        f"{phase_peak_yrfrac[pid]:.3f}", f"{r2[pid]:.3f}", f"{rms[pid]:.2f}",
                        f"{coh.reshape(-1)[pid]:.3f}" if coh is not None else "",
                        f"{rateA[pid]:.2f}", f"{rateB[pid]:.2f}",
                        f"{frac_down[pid]:.2f}" if np.isfinite(frac_down[pid]) else "",
                        cls_names[CLS[pid]]])

    # hotspot plots
    hs_rows = []
    if a.hotspots and os.path.exists(a.hotspots):
        import csv as csv2
        for row in csv2.DictReader(open(a.hotspots)):
            nm, la, lo = row["name"], float(row["lat"]), float(row["lon"])
            pub = row.get("published_rate_mm", "")
            j = int(round((lo - corner_lon) / post_lon - 0.5))
            i = int(round((la - corner_lat) / post_lat - 0.5))
            if not (0 <= i < length and 0 <= j < width):
                hs_rows.append((nm, la, lo, pub, None, None, "outside"))
                continue
            # 3x3 window mean
            i0, i1 = max(0, i - 1), min(length, i + 2)
            j0, j1 = max(0, j - 1), min(width, j + 2)
            block = cum[:, i0:i1, j0:j1].reshape(n_im, -1)
            wsel = ok_mask.reshape(sh)[i0:i1, j0:j1].reshape(-1)
            if wsel.sum() == 0:
                hs_rows.append((nm, la, lo, pub, None, None, "masked"))
                continue
            ts = np.nanmean(block[:, wsel], 1)
            pid = i * width + j
            bb, aa = b_rate[pid], amp_ann[pid]
            hs_rows.append((nm, la, lo, pub, bb, aa, "ok"))
            fig, ax = plt.subplots(figsize=(9, 4))
            ax.plot(dates_dt, ts, ".", ms=4, color="0.35", label="LiCSBAS cum. displacement")
            Gp = build_G(t, True)
            mfit = np.linalg.lstsq(Gp[np.isfinite(ts)], ts[np.isfinite(ts)], rcond=None)[0]
            tt = np.linspace(t[0], t[-1], 400)
            dd = [dates_dt[0] + (dates_dt[-1] - dates_dt[0]) * x / t[-1] for x in tt]
            ax.plot(dd, build_G(tt, True) @ mfit, "-", color="crimson", lw=1.6,
                    label=f"fit: {mfit[1]:.1f} mm/yr trend, {np.hypot(mfit[2],mfit[3]):.1f} mm annual amp")
            ax.plot(dd, mfit[0] + mfit[1] * tt, "--", color="navy", lw=1.2, label="secular trend")
            ax.set_title(f"{a.name} — {nm} ({la:.3f}N {lo:.3f}E)" + (f" | published ~{pub} mm/yr" if pub else ""))
            ax.set_ylabel("LOS displacement (mm)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
            fig.tight_layout()
            fig.savefig(out(f"hotspot_{nm.replace(' ', '_')}.png"), dpi=140)
            plt.close(fig)

    summary = {
        "name": a.name, "epochs": n_im, "span_yr": round(float(span), 2),
        "date_first": imdates[0], "date_last": imdates[-1],
        "usable_px": int(ok_mask.sum()), "total_px": int(length * width),
        "inelastic_px": int(inel.sum()), "elastic_px": int(elas.sum()),
        "indeterminate_px": int(valid.sum() - inel.sum() - elas.sum()),
        "rate_p2": float(np.nanpercentile(b_rate[valid], 2)),
        "rate_p50": float(np.nanpercentile(b_rate[valid], 50)),
        "rate_p98": float(np.nanpercentile(b_rate[valid], 98)),
        "amp_p50": float(np.nanpercentile(amp_ann[valid], 50)),
        "hotspots": [{"name": r[0], "lat": r[1], "lon": r[2], "published": r[3],
                      "fit_rate_mmyr": (None if r[4] is None else round(float(r[4]), 1)),
                      "fit_amp_mm": (None if r[5] is None else round(float(r[5]), 1)),
                      "status": r[6]} for r in hs_rows],
    }
    with open(out("summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))

if __name__ == "__main__":
    main()
