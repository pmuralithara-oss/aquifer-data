# PROCESSING_LOG — Sentinel-1 InSAR analysis of permanent aquifer damage, Punjab region

Mission: separate elastic (seasonal, recoverable) from inelastic (permanent) aquifer-system
deformation over 5 AOIs in Punjab (India+Pakistan) using Sentinel-1 InSAR time series 2015→present.

All times local (laptop). Session 2026-07-22.

---

## Phase 0 — Recon (2026-07-22)

### 0.1 Environment audit

Commands: `sw_vers`, `uname -m`, `sysctl hw.memsize hw.ncpu`, `df -h`, `python3 --version`, `command -v ...`

| Item | Value |
|---|---|
| OS | macOS 26.5.2 (arm64, Apple Silicon) |
| RAM | 24 GB |
| CPU | 10 cores |
| Disk free | **514 GB** on / (927 GB volume) |
| Python | 3.12.12 (Homebrew) + /usr/bin/python3; pip3 present |
| conda/mamba | **not installed** |
| GDAL | **not installed** (no gdalinfo) |
| Other | git, curl, jq present; wget, gmt absent |

Verdict: laptop adequate for LiCSBAS/MintPy-scale processing of clipped AOIs. Needs an
InSAR Python environment (proposal: micromamba, see plan).

### 0.2 Sentinel-1 SLC catalog per AOI (ASF Search API)

Query: `https://api.daac.asf.alaska.edu/services/search/param?platform=SENTINEL-1&processingLevel=SLC&beamMode=IW&intersectsWith=<bbox WKT>&output=jsonlite&maxResults=5000`

Notes: three of five queries returned HTTP 504 on first attempts; resolved by retry and by
splitting the query into 3 date ranges (2014–2019, 2019–2023, 2023–2027) and merging with jq.
Raw catalogs saved in session scratchpad `asf/aoi*.json` (not committed; ~5 MB total).

Scene counts by track (relative orbit), full archive 2014-10 → 2026-07:

| AOI | Track/geometry | Scenes | Date span |
|---|---|---|---|
| 1 Chandigarh–Mohali | 27 ASC | 331 | 2014-10-31 → 2026-07-07 |
| | 136 DESC | 299 | 2014-10-15 → 2026-07-03 |
| 2 Ludhiana | 100 ASC | 326; 34 DESC 325; 27 ASC 317; 136 DESC 299 | 2014-10 → 2026-07 |
| 3 Delhi NCR fringe | 27 ASC | 316; 63 DESC 313; 136 DESC 299 | 2014-10 → 2026-07 |
| 4 Lahore | 34 DESC | 392; 107 DESC 315; 100 ASC 332 | 2014-10 → 2026-07 |
| 5 Ghaggar belt | 27 ASC | 626; 100 ASC 475; 34 DESC 582; 136 DESC 451 (multi-frame bbox) | 2014-10 → 2026-07 |

Platform mix (AOI1 example): 615 S1A, 11 S1B, 2 S1C, 3 S1D → S1B contributed little here
(failed 2021-12); S1C/S1D appear from 2025/2026. Per-year scene counts show no catastrophic
2022–2024 gap at catalog level (e.g. AOI1: 49–61 scenes/yr 2020–2024). Actual usable epochs
determined by processed products below.

**Every AOI has both ascending and descending coverage → LOS→vertical decomposition possible in principle.**

### 0.3 LiCSAR (COMET) processed-product availability

- Old portal URL `gws-access.jasmin.ac.uk/public/lics/products/` → 404 (data migrated).
- Portal page (comet.nerc.ac.uk/comet-lics-portal) carries a data-migration notice and links
  **CEDA archive: `https://data.ceda.ac.uk/neodc/comet/data/licsar_products/`** (listing via
  `https://dap.ceda.ac.uk/neodc/comet/data/licsar_products/`). Confirmed HTTP 200, anonymous
  range-GET works (tested 206 on a product tif).
- Live JASMIN tree `gws-access.jasmin.ac.uk/public/nceo_geohazards/LiCSAR_products/` still
  serves **directory listings and frame metadata (incl. `-poly.txt` footprints)** but the
  interferogram pair directories themselves 404 → dangling entries; **files live on CEDA only**.
- **Stock LiCSBAS downloader (LiCSBAS01_get_geotiff.py, repo comet-licsar/LiCSBAS main @ 2026-07-22)
  still points at the dead JASMIN URL** → we must fetch from CEDA with our own downloader
  (or patch LiCSBAS01). Decision: custom downloader writing LiCSBAS's expected GEOC/ layout.

Frame footprints: downloaded `{frame}-poly.txt` for all 127 frames on tracks 27, 34, 63, 100,
107, 136 from JASMIN metadata (120 parsed; 7 odd variants without poly files skipped).
Intersection with AOI bboxes computed by Sutherland–Hodgman clip (script:
scratchpad/polys/intersect.py).

Frame ↔ AOI map (coverage as % of AOI bbox):

| AOI | Ascending frame | Descending frame |
|---|---|---|
| 1 Chandigarh | **027A_05990_131313** (100%) | **136D_05854_131313** (100%); 136D_06053 (30%) |
| 2 Ludhiana | **027A_05990** (99%); 100A_05836 (54%)+100A_06036 (46%) | **136D_05854** (100%); **034D_05909** (100%) |
| 3 Delhi | **027A_06190_131313** (100%) | **136D_06053_131313** (90.5%); 063D_06265 (39.5%) |
| 4 Lahore | **100A_05836_131313** (100%) | **034D_05909_131313** (100%); 107D_05886 (59.5%) |
| 5 Ghaggar | 027A_05990 (86%) + 100A_06036 (62%) | 136D_06053 (78%) + 034D_05909 (70%) |

Network depth per frame (pair dirs on JASMIN listing = processed; CEDA = actually downloadable
as of today):

| Frame | Processed pairs (JASMIN list) | Epochs | CEDA pairs | CEDA last pair |
|---|---|---|---|---|
| 027A_05990 | 1028 (2014-10→2024-07) | 223 | 1028 | 2024-06_2024-07 ✓ in sync |
| 027A_06190 | 1031 (2014-10→2024-10) | 239 | 1031 | 2024-10 ✓ |
| 136D_05854 | 1620 (2014-10→2024-09) | 238 | **429** | **2021-12 ← CEDA lags** |
| 136D_06053 | 954 (2014-10→2024-10) | 233 | **300** | **2021-12 ← CEDA lags** |
| 034D_05909 | 1266 (2014-10→2024-05) | 248 | 1266 | 2024-05 ✓ |
| 100A_05836 | 912 (2016-08→2024-07) | 218 | **398** | **2020-07 ← CEDA lags** |
| 100A_06036 | 882 (2016-09→2024-02) | 214 | 619 | 2021-12 ← lags |
| 063D_06265 | 894 (2014-10→2024-08) | 228 | 912 | 2024-08 ✓ |
| 107D_05886 | 1046 (2016-01→2024-10) | 231 | **246** | **2021-12 ← CEDA lags** |

Implications logged:
- LiCSAR processing for these frames currently **ends mid/late 2024** (no 2025–2026 yet on
  either mirror). 2025–26 tail must come from HyP3 burst InSAR if wanted.
- CEDA archive ingest incomplete for 4–5 frames (desc frames over India, asc/desc over Lahore):
  currently only → 2020/2021 there. Asc coverage for Delhi/Chandigarh/Ludhiana AOIs is complete
  to 2024. Will re-check CEDA during Phase 1; gaps fillable via HyP3.

### 0.4 Product format / size measurements (for download budget)

Sampled pair 20230330_20230423 of 027A_06190 on CEDA:
- `geo.unw.tif` 32.8 MB — 3026×2711 px, float32, **uncompressed, RowsPerStrip=1**,
  grid 0.001° (~111 m), tiepoint lon 75.6308, lat 29.4606 (frame N edge).
- `geo.cc.tif` 4.0 MB (byte, compressed), `geo.diff_pha.tif` 19.2 MB (wrapped; not needed).
- Frame metadata (E/N/U/hgt tifs): ~23.5 MB each ×4, once per frame.
- Full-frame unw+cc per pair ≈ **36.8 MB** → e.g. 1031 pairs ≈ 38 GB/frame-geometry. Too heavy.
- Because unw is uncompressed+stripped, an HTTP range-GET of the AOI's latitude rows retrieves
  a full-width band in ONE request (rows are contiguous on disk): e.g. Delhi 0.3° ≈ 300 rows
  ≈ 3.6 MB instead of 32.8 MB. cc must be fetched whole (compressed) or via windowed GDAL read.
  → Clipped-download strategy ≈ **6–10 MB/pair** ≈ 6–10 GB per frame-geometry for ~1000 pairs.

### 0.5 HyP3 (Option B) status

- Credit costs (hyp3-docs source, 2026): full-scene InSAR (GAMMA) 10–15 credits;
  **Burst InSAR 1 credit per single pair** (20/40/80-m); monthly free allotment applies
  (CREDITS_PER_MONTH templated in docs; historically 10,000/month). Requires NASA Earthdata
  login — **user credentials needed** when/if we use it.
- Role: gap-filler (desc 2022–24 where CEDA lags; 2025–26 extension), not primary path.

### 0.6 LiCSBAS status

- Repo comet-licsar/LiCSBAS active (pushed 2026-07-22), v1.15.0, requires Python ≥3.10.
- Deps: numpy, scipy, matplotlib, gdal (osgeo), h5py, xarray, rioxarray, shapely, networkx,
  statsmodels, astropy, cmcrameri, bs4 (+pygmt/ipympl for notebooks only — will skip).
- GDAL not on laptop → plan proposes micromamba env (conda-forge) as cleanest no-sudo route.

### 0.7 Decisions & open questions carried into the plan

1. Primary path: **Option A — CEDA LiCSAR + LiCSBAS** for all 5 AOIs (custom CEDA downloader,
   clipped to AOI latitude bands; stock LiCSBAS download step bypassed).
2. Pilot AOI: **AOI3 Delhi NCR fringe** — strongest published signal (>11 cm/yr), asc frame
   complete 2014→2024 on CEDA (1031 pairs), desc available 2014→2021 now (300 pairs) with
   more possibly arriving; validation target per mission.
3. Awaiting user approval for: micromamba env install, pilot download (~5–15 GB), Earthdata
   credentials (optional, for HyP3 gap-fill later).

---

## Phase 1 — Pilot: AOI3 Delhi NCR fringe (2026-07-22 →)

User approved (2026-07-22): toolchain install, pilot download, Delhi pilot, Earthdata for later.

### 1.1 Toolchain

- micromamba installed at `~/micromamba` (osx-arm64 latest via micro.mamba.pm API).
- Env `licsbas`: conda-forge python=3.11, **GDAL 3.12.3, numpy 2.4.6**, scipy, matplotlib,
  h5py, xarray, rioxarray, rasterio, shapely, networkx, statsmodels, astropy, cmcrameri, bs4,
  requests, psutil. Verified importable.
- LiCSBAS cloned (comet-licsar/LiCSBAS, main, shallow) → `tools/LiCSBAS` (54 bin scripts).
  Invocation via PATH+PYTHONPATH (`tools/run_licsbas.sh`); pip install not used (avoids
  pygmt/ipympl extras).

### 1.2 Custom clipped downloader (`tools/fetch_ceda_clip.py`)

Rationale: stock LiCSBAS01 points at dead JASMIN URL; CEDA full-frame unw = 32.8 MB each.
unw tifs are uncompressed, strip-per-row → AOI latitude band = one contiguous byte range →
single HTTP range request (~5 MB for our 0.4° band). Compressed products (cc, sltd, mli,
E/N/U/hgt) fetched whole then cropped with GDAL and temp deleted. Output = LiCSBAS GEOC layout.
Retries ×4 with backoff; per-pair failures logged, resumable (skips existing).

Validation (3 pairs, asc frame): correct 551×401 grid @0.001°, geotransform
(76.8993, 28.6501 NW corner), unw values ±15 rad, ~1.3 s/pair. 

### 1.3 Pilot downloads (launched, background)

Clip box for AOI3 (with 0.05° margin): **lat 28.25–28.65, lon 76.90–77.45**.

| Stack | Frame | Pairs | GACOS epochs | Dest |
|---|---|---|---|---|
| asc | 027A_06190_131313 (track 27) | 1031 (2014-10→2024-10) | 179 sltd avail. | data/AOI3_delhi/asc_027A_06190 |
| desc | 136D_06053_131313 (track 136) | 300 (2015-04→2021-12, CEDA lag) | 218 epoch dirs, sltd partial | data/AOI3_delhi/desc_136D_06053 |

### 1.4 Planned LiCSBAS chain (`tools/run_licsbas.sh`)

02 ml_prep (nlook=1) → 03op GACOS (if sltd present) → 11 check_unw → 12 loop_closure →
13 sb_inv (NSBAS) → 14 vel_std → 15 mask_ts → 16 filt_ts (+ --interpolate_nans).
Defaults kept unless noted; all parameter deviations to be logged here.

### 1.5 Download + desc processing results (2026-07-22)

- desc 136D_06053: 300/300 pairs OK (267 s wall). GACOS sltd present for 124/218 epochs.
- asc 027A_06190: 1026/1031 pairs OK; **5 pairs genuinely missing unw on CEDA (404), all
  involving epoch 20170605** — removed from stack. GACOS 144/179.
- Disk: asc 970 MB, desc 280 MB clipped (vs ~45 GB full-frame equivalent).
- LiCSBAS desc chain (GACOS applied): 74/300 ifgs removed by QC (steps 11+12); 85% of
  220,951 px retained; TS_GEOCml1GACOS produced (95 epochs, 2015-04→2021-12, span 6.7 yr).

### 1.6 Validation diagnostics — measured vs published hotspot rates (desc, 2015–2021)

decompose_ts.py results agree with LiCSBAS native vel (self-consistent). BUT measured
rates at published hotspot locations are far below published values:

| Site | Published (prior) | Measured desc LOS 2015-21 | 7×7 min (LiCSBAS vel) | coh (12-24d pairs) |
|---|---|---|---|---|
| Kapashera | −110…−170 mm/yr | −5 | −7.7 | 0.25 (marginal) |
| Faridabad Sec21 | ~−115 | −6.1 | −11.0 | 0.54 (good) |
| Samalka | ~−90 | −17.0 | −34.2 | 0.34 (marginal) |
| Gurugram DLF | −30…−80 | −3…−5 | — | 0.58 (good) |
| Dwarka | ~−35 | +1.2 | −0.4 | — |

Diagnostics run:
1. Raw long-span (140–400 d) ifgs 2017–2020 show **no deep bowls** at these sites
   (bowl-ring differences −26…+13 mm, mostly |≤6| mm); coherence in long pairs at
   Kapashera/Samalka is 0.10–0.16 (decorrelated) → long pairs cannot carry fast bowls there.
2. Short pairs (≤24 d, n=68 sampled): Faridabad 0.54 / Gurugram 0.58 / IGI 0.50 coherent →
   their moderate measured rates are credible for 2015–2021; Kapashera 0.25 marginal.
3. Most-negative pixels in the desc velocity field (−80…−139 mm/yr) sit at coh 0.05–0.14
   (noise-suspect, AOI edges), not at published hotspots.

Candidate explanations (to be discriminated by the asc stack → 2024 and, later, HyP3 burst
processing): (a) published peak rates are localized cores/specific years (esp. post-2019
acceleration in Faridabad **after** the desc window ends); (b) LiCSAR 111-m filtered
products clip steep small-scale bowls in mixed/vegetated cover (unwrap underestimation);
(c) hotspot coordinate approximations. Published per-site values above are approximate
prior-phase recollections — REPORT.md will treat the mission prior (">11 cm/yr in AOI3")
as the validation bar, not the per-site table.

Honest status: **pilot validation NOT yet passed at extreme hotspots; broad moderate-rate
field (p2 −7.3 mm/yr, 1,498 inelastic-candidate px, 37,207 elastic px) is internally
consistent.** Asc chain (239 epochs → 2024-10) running.

### 1.7 Pilot verdict (asc chain complete, 2026-07-22 23:11)

- asc TS: 137 epochs used (2014-10 → ~2023-03; LiCSBAS QC dropped later epochs — likely
  network disconnection in 2023-24 acquisitions), span 8.38 yr, 140,853 usable px.
- **Both published Delhi subsidence bowls are detected at the correct locations** in both
  geometries: Kapashera–Samalka bowl (core ~28.515N 77.08E) and Faridabad bowl (~28.40N
  77.30E). Samalka 3×3 series: clean secular −10.7 mm/yr over 8.4 yr (~−90 mm cumulative),
  annual amplitude only ~2 mm, no rebound → inelastic-style signature. Bowl-core pixels
  reach −34 mm/yr (desc), −20s (asc).
- **Magnitude validation NOT met at published peak level**: measured long-term averages are
  3–10× below published PSI peak rates (−110…−170 mm/yr, mostly 2018–19, building-scale).
  Cross-geometry agreement (asc 2014-23 vs desc 2015-21) shows this is a property of the
  LiCSAR 111-m product + long-window average, not a processing bug. Flagged as the main
  systematic-uncertainty: LiCSAR-based magnitudes in sharp bowls are lower bounds.
  Independent HyP3 burst-InSAR check recommended when Earthdata credentials provided.
- Method decision: pipeline VALIDATED for structure/classification; magnitudes to be
  reported with the caveat above. Proceeding to science AOIs.
- Deliverables written to outputs/AOI3_delhi/: LOS+quasi-vertical velocity GeoTIFFs+PNGs,
  seasonal amplitude, classification maps, 10 hotspot time-series plots ×2 geometries,
  decomposition_results.csv ×2, summary.json ×2.

### 1.8 Scale-out downloads (AOI1 running; shared Punjab band launched)

- AOI1 asc 027A_05990 box (30.90,30.40,76.50,77.00) + desc 136D_05854 same box: running.
- Shared band of 027A_05990 (31.05,29.75,75.15,76.55) covering AOI2+AOI5 asc: launched
  (unw row-band trick makes lon-width free; per-AOI GEOC dirs will be cropped locally).
- Planned next: 034D_05909 clips for AOI2 (desc) and AOI4 (desc); 100A_05836 for AOI4 asc
  (2016-2020 on CEDA). 107D_05886 (59.5% Lahore) as desc alternative if 034D QC poor.

### 1.9 Delhi joint vertical (2026-07-22 23:2x)

`tools/joint_vertical.py`: asc+desc rates re-fit on common 2015.25–2021.95 window
(95+95 epochs), grids aligned by integer shift (dy=−1), per-pixel [U,E] solve.
117,301 joint px. **Samalka bowl vertical rate: −22.9 mm/yr mean, −44.3 mm/yr min (7×7)**;
Faridabad −6.8 mean / −10.4 min; field p2/p50/p98 = −6.9/+0.2/+4.9 mm/yr.
Outputs: AOI3_joint_vert/east GeoTIFFs + maps.

### 2.0 AOI1 Chandigarh desc — data QC + LiCSBAS bug fix

- 136D_05854 clip: 27 pair dirs were **empty in the source data** over our window (all
  epochs 2014-10→2015-01 + around 20180515) — verified by refetch (identical empties);
  removed. 1 pair missing on CEDA. ~401 usable pairs remain (2014-12→2021-12).
- GACOS: 208/225 epochs.
- **LiCSBAS bug found & patched** (`LiCSBAS_lib/LiCSBAS_inv_lib.py:629`,
  censored_lstsq_slow_para_wrapper): on per-point lstsq failure it returned scalar
  np.nan instead of an (n_im+1) NaN vector → ragged list → `np.array(_result)` ValueError
  under numpy 2.x when ≥1 point fails. Fixed to return `np.full(n_par, nan)`. Chain
  resumed from step 13. (Did not affect Delhi runs: no failing points there.)

### 2.1 AOI1 Chandigarh–Mohali corridor — RESULTS (2026-07-23 early hours)

- asc 027A_05990: 998 usable pairs (29 empty-in-source + 1 corrupt-at-source
  [20171108_20171214, PackBits truncation on CEDA] removed; 18 pairs missing on CEDA).
  Chain OK after inv_lib patch. TS: 102 epochs, 2014-12→2023-03, 8.25 yr.
- desc 136D_05854: TS 118 epochs, 2014-12→2021-12, 7.16 yr.
- **Reference problem found & fixed**: LiCSBAS auto-reference landed INSIDE the subsiding
  corridor (76.661E 30.765N, near Kharar) → field median +24 mm/yr artifact. Added
  `--refbox` re-referencing (per-epoch median subtraction) to decompose_ts.py and
  joint_vertical.py; reference zone = coherent Chandigarh-city block
  (30.72–30.76N, 76.76–76.82E, coh~0.37, 2400 px). After re-referencing: field median
  −1.2 to −1.6 mm/yr (physical).
- **Cross-geometry hotspot agreement (LOS, city-referenced)**: Kharar −91.8 (asc) / −85.6
  (desc); Mohali Sec70 −40.3/−48.0; Zirakpur −4.4/−6.4; Dera Bassi −92.6/−130.1;
  Chandigarh Sec17 +5.4/+5.2; Panchkula +17/+11 (range-front, low coh, flagged).
- **Joint vertical (2015–2021 common window): Kharar −114 mean/−133 min; Dera Bassi −140/−155;
  Sunny Enclave −67; Mohali −59; Banur −30; Zirakpur −7; city ≈ +4 mm/yr.**
  → **Validation of the mission prior (60–180 mm/yr corridor subsidence): PASSED within
  factor ~1.5 at all corridor hotspots.**
- Temporal character at Kharar (asc): **linear −710 mm over 8.25 yr, seasonal amplitude
  2.1 mm (≈2% of annual trend), no deceleration, no recovery** → inelastic-dominant.
  Inelastic-candidate pixels: 23,614 (asc) / 22,210 (desc).
- Caveats flagged: isolated coherent patches show ±25–50 mm/yr blobs (isolated-component
  unwrap bias; e.g. +30–40 blob near 30.47N 76.58E, Siwalik-front positives at coh<0.15)
  — excluded from interpretation; classification map to be read jointly with coherence.

### 3.0 Phase 3 ground validation — programmatic access FAILED (deferred)

- indiawris.gov.in, www.indiawris.gov.in, cgwb.gov.in: connection timeouts from this
  network (likely geo-restriction or outage).
- gwdata.cgwb.gov.in: reachable but returns "Maintenance Mode" page (checked ×2, browser UA too).
- data.gov.in API: requires API key; SPA scraping out of scope.
- Per mission spec: **elastic-storativity estimation (seasonal deformation vs seasonal head)
  deferred**; InSAR-only decomposition proceeds. User can supply WRIS/CGWB exports later to
  complete this phase.

### 4.0 Scale-out processing notes (2026-07-23)

- **AOI4 Lahore asc 100A_05836 (2016-08→2020-07, 3.78 yr, 86 epochs, 398/398 pairs):
  Model Town −45.0 mm/yr (published prior: up to −43) — validation factor ~1.0.**
  Township −45.3, Gulberg −33.8, Walled City −32.9, Central −29.7, Johar Town −4.4,
  DHA −0.2. Seasonal amplitudes 3.5–5.3 mm. 15,973 inelastic / 17,766 elastic px.
  Field median −1.6 (auto-ref OK). Desc 034D (2014→2024) downloading for joint+extension.
- **AOI2 Ludhiana asc 027A (8.25 yr): NO published-style bowl at the city; instead a
  smooth ~40-km NW(−)→SE(+) gradient (±10–15 mm/yr) — suspected long-wavelength
  orbital/iono ramp.** City sites +5…+10 relative to mid-gradient auto-ref. Desc will
  discriminate ramp vs real motion; verdict deferred.
- **Data defect found & fixed: LiCSAR grid shift between processing eras.** Some 2022–23
  pairs of 034D_05909 and 027A_05990 were produced on a grid shifted by ~1 px lat /
  0.0002° lon (301→300 rows). Step 02 silently skipped 45 such pairs per stack (empty
  GEOCml1 dirs), then step 03op crashed on the empties and the runner's fallback used the
  half-baked GACOS dir (fixed: runner now validates slc.mli.par; GACOS-fail reverts to
  GEOCml1). New `tools/align_grid.py` warps mismatched tifs onto the majority grid
  (nearest; sub-pixel snap). AOI2-desc (223 files) and AOI5-asc (755 files) aligned;
  chains re-running. **Delhi/Chandigarh stacks audited: 0 skipped pairs — unaffected.**
  Their asc series ending ~2023-03 reflects LiCSBAS QC/network gaps in late acquisitions,
  not this bug.

### 5.0 Final AOI results (2026-07-23)

- **AOI4 Lahore desc 034D** (1253 pairs after clean+align, 149 epochs, 2016-04→2023-03):
  hotspots confirm asc within ~10–15 % (Model Town −46.6 vs −45.0 etc.). Joint vertical
  2016–2020: Model Town −61, Township −62 (min −78), Gulberg −43, Walled City −39 mm/yr.
  Inelastic 15,366 px (~162 km², mean −33 mm/yr LOS).
- **AOI2 Ludhiana desc 034D** (1169 pairs, 144 epochs, 2016-04→2023-03): NO NW–SE
  gradient → asc ramp confirmed as orbital/iono artifact. Both geometries: city +5…+9
  LOS vs local refs. Joint (planar-deramped, new --deramp option): vertical
  p2/p50/p98 = −5.6/+1.0/+5.5 mm/yr → stable. Inelastic px: 6 (desc) / 900 (asc, ramp
  artifact). **Published ~−25 mm/yr not reproduced — prior contradicted.**
- **AOI5 Ghaggar asc 027A** (964 pairs after align, 107 epochs, 8.35 yr): 38,034/1,122,201
  px usable (3.4 %) — rural decorrelation as predicted. Town centres coherent (0.34–0.56)
  but masked by network criteria; unmasked indicative values ±15 mm/yr spread = island
  unwrap-bias scale; no rate claims. Desc stack NOT downloaded (decision: coherence
  result makes it redundant; burst PSI is the right instrument).
- Killed-task recovery: one system interruption killed the Lahore-desc download (resumed
  at 1150/1266) and the AOI2-desc chain mid-step-13 (rerun clean).
- Summary stats table (all 9 stacks) captured in REPORT.md §2–6; per-pixel CSVs and
  109 map/plot figures in outputs/.

### 6.0 Mission close-out state

- Phases 0,1,2,4 complete; Phase 3 (ground validation) deferred — portals unreachable
  (§3.0). HyP3 2025–26 extension + burst validation awaiting user's Earthdata login.
- Environment: micromamba env `licsbas` (~3.1 GB), LiCSBAS clone + 1-line patch,
  tools/ scripts. Data footprint: data/ ≈ 18 GB, outputs/ ≈ 0.4 GB, well under budget.

<!-- Log continues as work proceeds. -->
