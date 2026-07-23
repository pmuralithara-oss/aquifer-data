#!/usr/bin/env python3
"""
Clipped-window downloader for COMET-LiCSAR products on CEDA.

Why: full-frame geo.unw.tif files are ~33 MB but uncompressed with RowsPerStrip=1,
so the rows covering an AOI latitude band are one contiguous byte range -> we fetch
them with a single HTTP range request (~5 MB) and write a windowed GeoTIFF locally.
Compressed files (cc, sltd, mli, E/N/U/hgt) are fetched whole to a temp file, cropped
with GDAL, and the temp deleted.

Output follows the LiCSBAS GEOC/ layout so LiCSBAS can run from step 02:
  OUTDIR/GEOC/{pair}/{pair}.geo.unw.tif , .geo.cc.tif
  OUTDIR/GEOC/{frame}.geo.E.tif .N .U .hgt .mli.tif, baselines, metadata.txt
  OUTDIR/GACOS/{date}.sltd.geo.tif
"""
import argparse, io, os, shutil, struct, sys, tempfile, threading, time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import requests
from osgeo import gdal, osr

gdal.UseExceptions()
CEDA = "https://dap.ceda.ac.uk/neodc/comet/data/licsar_products"
TIMEOUT = 120
HEAD_BYTES = 262143  # 256 KB: covers IFD + strip offset arrays for these files

_tls = threading.local()

def sess():
    if not hasattr(_tls, "s"):
        s = requests.Session()
        s.headers["User-Agent"] = "aquifer-insar-research/0.1 (LiCSBAS prep; range reads to minimise load)"
        _tls.s = s
    return _tls.s

def get_range(url, start, end, tries=4):
    for a in range(tries):
        try:
            r = sess().get(url, headers={"Range": f"bytes={start}-{end}"}, timeout=TIMEOUT)
            if r.status_code in (200, 206):
                return r.content
            if r.status_code == 404:
                raise FileNotFoundError(url)
        except FileNotFoundError:
            raise
        except Exception:
            pass
        time.sleep(2 * (a + 1))
    raise RuntimeError(f"range fetch failed: {url}")

def get_whole(url, dst, tries=4):
    for a in range(tries):
        try:
            with sess().get(url, stream=True, timeout=TIMEOUT) as r:
                if r.status_code == 404:
                    raise FileNotFoundError(url)
                r.raise_for_status()
                with open(dst, "wb") as fh:
                    shutil.copyfileobj(r.raw, fh)
            return
        except FileNotFoundError:
            raise
        except Exception:
            time.sleep(2 * (a + 1))
    raise RuntimeError(f"download failed: {url}")

TYPE_SIZE = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 7: 1, 11: 4, 12: 8}

def parse_tiff_header(url):
    """Parse enough of a classic little-endian TIFF to plan a strip-band range read."""
    d = get_range(url, 0, HEAD_BYTES)
    if d[:2] != b"II" or struct.unpack("<H", d[2:4])[0] != 42:
        return None  # BigTIFF or big-endian: caller falls back to whole-file
    off = struct.unpack("<I", d[4:8])[0]
    if off + 2 > len(d):
        return None
    n = struct.unpack("<H", d[off:off + 2])[0]
    tags = {}
    need_more = 0
    for i in range(n):
        e = d[off + 2 + i * 12: off + 2 + (i + 1) * 12]
        if len(e) < 12:
            return None
        tag, typ, cnt = struct.unpack("<HHI", e[:8])
        size = TYPE_SIZE.get(typ, 1) * cnt
        if size <= 4:
            if typ == 3:
                val = list(struct.unpack(f"<{cnt}H", e[8:8 + 2 * cnt]))
            elif typ == 4:
                val = list(struct.unpack(f"<{cnt}I", e[8:8 + 4 * cnt]))
            else:
                val = [struct.unpack("<I", e[8:12])[0]]
            tags[tag] = ("inline", typ, cnt, val)
        else:
            ptr = struct.unpack("<I", e[8:12])[0]
            tags[tag] = ("ptr", typ, cnt, ptr)
            need_more = max(need_more, ptr + size)
    if need_more > len(d):
        d = d + get_range(url, len(d), need_more - 1)

    def read_tag(tag):
        if tag not in tags:
            return None
        kind, typ, cnt, v = tags[tag]
        if kind == "inline":
            return v
        ptr = v
        size = TYPE_SIZE[typ] * cnt
        raw = d[ptr:ptr + size]
        fmt = {3: "H", 4: "I", 12: "d", 2: "s", 11: "f"}[typ]
        if typ == 2:
            return raw.rstrip(b"\x00").decode("ascii", "replace")
        return list(struct.unpack(f"<{cnt}{fmt}", raw))

    info = {
        "width": read_tag(256)[0], "height": read_tag(257)[0],
        "bps": read_tag(258)[0], "compression": read_tag(259)[0],
        "sf": (read_tag(339) or [1])[0],
        "rps": (read_tag(278) or [None])[0],
        "strip_offsets": read_tag(273), "strip_counts": read_tag(279),
        "scale": read_tag(33550), "tiepoint": read_tag(33922),
        "nodata": read_tag(42113),
    }
    return info

def np_dtype(info):
    bps, sf = info["bps"], info["sf"]
    if sf == 3 and bps == 32: return np.float32
    if sf == 1 and bps == 8: return np.uint8
    if sf == 2 and bps == 32: return np.int32
    if sf == 1 and bps == 16: return np.uint16
    raise ValueError(f"unhandled sample format {sf}/{bps}")

GDT = {np.dtype(np.float32): gdal.GDT_Float32, np.dtype(np.uint8): gdal.GDT_Byte,
       np.dtype(np.int32): gdal.GDT_Int32, np.dtype(np.uint16): gdal.GDT_UInt16}

def write_tif(path, arr, gt, nodata=None):
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, arr.shape[1], arr.shape[0], 1, GDT[arr.dtype],
                    options=["COMPRESS=DEFLATE", "PREDICTOR=1", "TILED=YES"])
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference(); srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1)
    if nodata is not None:
        try: b.SetNoDataValue(float(nodata))
        except Exception: pass
    b.WriteArray(arr)
    ds.FlushCache(); ds = None

def window_from_info(info, lat_n, lat_s, lon_w, lon_e):
    lon0, lat0 = info["tiepoint"][3], info["tiepoint"][4]
    dx, dy = info["scale"][0], info["scale"][1]
    x0 = max(0, int((lon_w - lon0) / dx)); x1 = min(info["width"], int(np.ceil((lon_e - lon0) / dx)))
    y0 = max(0, int((lat0 - lat_n) / dy)); y1 = min(info["height"], int(np.ceil((lat0 - lat_s) / dy)))
    if x1 <= x0 or y1 <= y0:
        raise ValueError("window outside raster")
    gt = (lon0 + x0 * dx, dx, 0.0, lat0 - y0 * dy, 0.0, -dy)
    return x0, x1, y0, y1, gt

def fetch_band_uncompressed(url, info, win, out, nodata_default=None):
    """Single range request for rows y0:y1 of an uncompressed strip-per-row tiff."""
    x0, x1, y0, y1, gt = win
    so, sc = info["strip_offsets"], info["strip_counts"]
    rowbytes = info["width"] * (info["bps"] // 8)
    for i in range(y0, y1 - 1):
        if so[i + 1] != so[i] + rowbytes:
            return False  # not contiguous -> caller falls back
    raw = get_range(url, so[y0], so[y1 - 1] + rowbytes - 1)
    arr = np.frombuffer(raw, dtype=np_dtype(info)).reshape(y1 - y0, info["width"])[:, x0:x1].copy()
    nd = info["nodata"] if info["nodata"] is not None else nodata_default
    write_tif(out, arr, gt, nodata=float(nd) if nd not in (None, "") else None)
    return True

def fetch_whole_and_crop(url, out, lat_n, lat_s, lon_w, lon_e, nodata_default=None):
    fd, tmp = tempfile.mkstemp(suffix=".tif", dir=os.path.dirname(out))
    os.close(fd)
    try:
        get_whole(url, tmp)
        ds = gdal.Open(tmp)
        gt0 = ds.GetGeoTransform()
        W, H = ds.RasterXSize, ds.RasterYSize
        x0 = max(0, int((lon_w - gt0[0]) / gt0[1])); x1 = min(W, int(np.ceil((lon_e - gt0[0]) / gt0[1])))
        y0 = max(0, int((gt0[3] - lat_n) / -gt0[5])); y1 = min(H, int(np.ceil((gt0[3] - lat_s) / -gt0[5])))
        if x1 <= x0 or y1 <= y0:
            raise ValueError(f"window outside raster for {url}")
        arr = ds.GetRasterBand(1).ReadAsArray(x0, y0, x1 - x0, y1 - y0)
        nd = ds.GetRasterBand(1).GetNoDataValue()
        gt = (gt0[0] + x0 * gt0[1], gt0[1], 0.0, gt0[3] + y0 * gt0[5], 0.0, gt0[5])
        ds = None
        write_tif(out, arr, gt, nodata=nd if nd is not None else nodata_default)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

def do_pair(track, frame, pair, box, geoc):
    lat_n, lat_s, lon_w, lon_e = box
    pdir = os.path.join(geoc, pair)
    os.makedirs(pdir, exist_ok=True)
    base = f"{CEDA}/{track}/{frame}/{pair}/{pair}"
    results = []
    for kind, nodata_default in (("geo.unw.tif", 0.0), ("geo.cc.tif", 0.0)):
        out = os.path.join(pdir, f"{pair}.{kind}")
        if os.path.exists(out) and os.path.getsize(out) > 2000:
            results.append("cached"); continue
        url = f"{base}.{kind}"
        try:
            done = False
            info = parse_tiff_header(url)
            if info and info["compression"] == 1 and info["rps"] == 1 and info["strip_offsets"]:
                win = window_from_info(info, lat_n, lat_s, lon_w, lon_e)
                done = fetch_band_uncompressed(url, info, win, out, nodata_default)
            if not done:
                fetch_whole_and_crop(url, out, lat_n, lat_s, lon_w, lon_e, nodata_default)
            results.append("ok")
        except FileNotFoundError:
            results.append(f"MISSING:{kind}")
        except Exception as e:
            results.append(f"FAIL:{kind}:{e}")
    return pair, results

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", required=True)
    ap.add_argument("--frame", required=True)
    ap.add_argument("--pairs", required=True, help="file with one pair per line")
    ap.add_argument("--box", required=True, help="latN,latS,lonW,lonE")
    ap.add_argument("--out", required=True, help="output dir (GEOC/, GACOS/ created inside)")
    ap.add_argument("--gacos", action="store_true", help="also fetch epoch sltd files")
    ap.add_argument("--epochs", help="file with one epoch date per line (for --gacos)")
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    lat_n, lat_s, lon_w, lon_e = [float(x) for x in a.box.split(",")]
    box = (lat_n, lat_s, lon_w, lon_e)
    geoc = os.path.join(a.out, "GEOC")
    os.makedirs(geoc, exist_ok=True)

    # frame-level metadata
    meta_dir = f"{CEDA}/{a.track}/{a.frame}/metadata"
    for name in ("metadata.txt", "baselines"):
        dst = os.path.join(geoc, name)
        if not os.path.exists(dst):
            try: get_whole(f"{meta_dir}/{name}", dst)
            except Exception as e: print(f"WARN metadata {name}: {e}")
    master = None
    mt = os.path.join(geoc, "metadata.txt")
    if os.path.exists(mt):
        for line in open(mt):
            if line.startswith("master="):
                master = line.strip().split("=")[1]
    for suff in ("geo.E.tif", "geo.N.tif", "geo.U.tif", "geo.hgt.tif"):
        out = os.path.join(geoc, f"{a.frame}.{suff}")
        if not os.path.exists(out):
            try:
                fetch_whole_and_crop(f"{meta_dir}/{a.frame}.{suff}", out, *box)
                print(f"meta {suff} done")
            except Exception as e:
                print(f"WARN {suff}: {e}")
    if master:
        out = os.path.join(geoc, f"{a.frame}.geo.mli.tif")
        if not os.path.exists(out):
            try:
                fetch_whole_and_crop(f"{CEDA}/{a.track}/{a.frame}/epochs/{master}/{master}.geo.mli.tif", out, *box)
                print("master mli done")
            except Exception as e:
                print(f"WARN master mli: {e}")

    pairs = [l.strip() for l in open(a.pairs) if l.strip()]
    if a.limit: pairs = pairs[:a.limit]
    print(f"{len(pairs)} pairs -> {geoc}  box={box}")
    t0 = time.time()
    ok = cached = 0
    fails = []
    with ThreadPoolExecutor(a.nproc) as ex:
        futs = {ex.submit(do_pair, a.track, a.frame, p, box, geoc): p for p in pairs}
        for i, fut in enumerate(as_completed(futs), 1):
            pair, res = fut.result()
            if all(r in ("ok", "cached") for r in res):
                ok += 1
                if all(r == "cached" for r in res): cached += 1
            else:
                fails.append((pair, res))
            if i % 50 == 0 or i == len(pairs):
                el = time.time() - t0
                print(f"[{i}/{len(pairs)}] ok={ok} fail={len(fails)} {el:.0f}s", flush=True)
    for p, r in fails:
        print("FAILED", p, r)

    if a.gacos and a.epochs:
        gac = os.path.join(a.out, "GACOS")
        os.makedirs(gac, exist_ok=True)
        eps = [l.strip() for l in open(a.epochs) if l.strip()]
        print(f"GACOS: {len(eps)} epochs")
        def do_ep(d):
            out = os.path.join(gac, f"{d}.sltd.geo.tif")
            if os.path.exists(out) and os.path.getsize(out) > 2000:
                return d, "cached"
            url = f"{CEDA}/{a.track}/{a.frame}/epochs/{d}/{d}.sltd.geo.tif"
            try:
                info = parse_tiff_header(url)
                done = False
                if info and info["compression"] == 1 and info["rps"] == 1 and info["strip_offsets"]:
                    win = window_from_info(info, *box)
                    done = fetch_band_uncompressed(url, info, win, out)
                if not done:
                    fetch_whole_and_crop(url, out, *box)
                return d, "ok"
            except FileNotFoundError:
                return d, "missing"
            except Exception as e:
                return d, f"fail:{e}"
        stats = {}
        with ThreadPoolExecutor(a.nproc) as ex:
            for d, st in ex.map(do_ep, eps):
                stats[st.split(":")[0]] = stats.get(st.split(":")[0], 0) + 1
        print("GACOS result:", stats)
    print(f"total wall {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
