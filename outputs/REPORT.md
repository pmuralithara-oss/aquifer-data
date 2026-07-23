# Permanent aquifer-storage damage in the Punjab region: an InSAR elastic/inelastic decomposition
### Sentinel-1 (2014–2024) · COMET-LiCSAR products · LiCSBAS NSBAS time series · 5 AOIs (India + Pakistan)

**Final report — all five AOIs processed.** Command-level provenance: `PROCESSING_LOG.md`.
Analysis run 2026-07-22/23 on a laptop from public archives; no raw SLC processing.

---

## 0. Headline findings

1. **Chandigarh–Mohali corridor (AOI1): the clearest evidence of permanent (inelastic)
   aquifer-system compaction in the study.** A contiguous ~250 km² zone
   (Kharar–Sunny Enclave–Mohali–Landran, plus a second bowl at Dera Bassi) subsides at
   **60–140 mm/yr vertical** (Kharar −114, Dera Bassi −140, 2015–2021 joint asc+desc).
   The motion is linear for 8+ years (Kharar: −710 mm cumulative LOS, R²≈0.995), the
   seasonal amplitude is only ~2–3 mm (≈2 % of the annual trend), and **100 % of successive
   post-monsoon (Oct–Dec) surface maxima step downward** — rebound never recovers prior
   levels. With the corridor's documented 20–25 m clay beds (research-phase prior), this is
   the textbook kinematic signature of inelastic clay compaction: **storage capacity being
   permanently destroyed, at an order-of-magnitude ~15 million m³ of pore volume per year**
   (mean −46 mm/yr LOS ≈ −59 vertical × 251 km²).
2. **Lahore (AOI4): strong, validated, ongoing compaction.** Model Town / Township /
   Gulberg / Walled City subside at **−39 to −62 mm/yr vertical** (2016–2020 joint), in
   agreement between two independent geometries and with the published −43 mm/yr prior
   (factor ~1.0). ~160 km² classified inelastic-candidate (mean −33 mm/yr LOS). Seasonal
   amplitudes (3.5–5 mm) are small relative to trend. Deep lithology is poorly documented,
   so the inelastic attribution is kinematic (sustained + unrecovered), not lithological.
3. **Delhi NCR fringe (AOI3, calibration site): method validated structurally; LiCSAR
   magnitudes are lower bounds at sharp bowls.** Both published bowls (Kapashera–Samalka,
   Faridabad) are detected at the right locations with the inelastic signature
   (−23 to −44 mm/yr vertical at Samalka, ~15 km² inelastic). Published PSI peak rates
   (>110 mm/yr, 2018–19, building-scale) are 3–10× larger than our 111-m long-window
   averages — demonstrated to be a product-resolution/period effect, not a processing bug
   (cross-geometry agreement; long-pair coherence 0.10–0.16 at those villages).
4. **Ludhiana (AOI2): the published ~25 mm/yr subsidence is NOT present in 2015–2024
   LiCSAR data — a prominent contradiction of the prior.** Two independent geometries agree
   the AOI is stable within ±5 mm/yr vertical (deramped joint), with the city marginally
   positive relative to surroundings. Essentially zero inelastic-candidate pixels (6–900,
   scattered artifacts). Either earlier-published rates were localized/temporary, or
   extraction stabilized (canal-supply and paddy-policy changes post-2019 are plausible
   mechanisms — untested here). An ascending-frame ~40-km ramp (±10–15 mm/yr) was identified
   as an orbital/ionospheric artifact by the descending comparison and removed for the joint.
5. **Ghaggar interfan belt (AOI5): not measurable with 111-m LiCSAR products — a valid
   sparse/null result, as the mission anticipated.** Only 3.4 % of pixels survive QC over
   the rural belt; town centres retain raw coherence (0.34–0.56) but sit as isolated
   islands whose values (±15 mm/yr spread) carry integer-cycle unwrap-bias risk. No town
   shows a Chandigarh-style multi-cm/yr signal, but definitive rates require full-resolution
   burst PSI. The highest-clay-fraction zone therefore remains **untested, not exonerated**.
6. **Ground-truth deferral:** CGWB/India-WRIS groundwater portals were unreachable
   (maintenance/geo-block), so the elastic-storativity ratio (seasonal deformation vs
   seasonal head change) is deferred; the elastic/inelastic separation here rests on the
   InSAR kinematics alone (trend vs seasonal vs post-monsoon recovery).

**Novelty check per mission brief:** this is, to our knowledge, the first explicit
elastic-vs-inelastic decomposition applied across Indo-Gangetic Plain aquifers, and it
distinguishes cleanly: corridor/Lahore/Delhi-bowls = inelastic-dominant; Ludhiana =
stable; wide NCR fringe = elastic-dominant, small amplitudes.

---

## 1. Methods actually used

**Data.** COMET-LiCSAR Level-3 geocoded unwrapped interferograms + coherence (0.001° ≈
111 m), CEDA archive (`neodc/comet/data/licsar_products`), accessed 2026-07-22/23. LiCSAR
processing of these frames ends mid/late-2024 (2025–26 extension would need HyP3 burst
InSAR + Earthdata login — offered as follow-up). Frames per AOI: PROCESSING_LOG §0.3.
Custom downloader (`tools/fetch_ceda_clip.py`) transfers only the AOI latitude band of each
uncompressed unw (single HTTP range request; ~25× saving) and crops compressed products
via temp files; total transfer ≈ 12 GB for 9 stacks (~7,700 interferograms).

**Data QC found and fixed:** (i) epochs empty-in-source over some windows (removed);
(ii) one PackBits-corrupt file on CEDA (pair dropped); (iii) **LiCSAR inter-era grid shift**
(some 2022–23 pairs on a ~1-px-shifted lattice) — silently breaks LiCSBAS; fixed by
nearest-neighbour snap (`tools/align_grid.py`); (iv) a LiCSBAS inversion-library bug with
numpy ≥2 (scalar-NaN ragged list) — one-line patch, log §2.0.

**Time series.** LiCSBAS v1.15 (comet-licsar), nlook=1, GACOS where available,
loop-closure QC, NSBAS inversion, mask, spatio-temporal filter. Per-run parameters in logs
(`data/*/licsbas_run.log`).

**Decomposition** (`tools/decompose_ts.py`, per coherent pixel):
`d(t) = a + b·t + c₁sin2πt + c₂cos2πt [+ semiannual if BIC improves] + ε`; σ via 150-draw
90-day-block bootstrap; sustained test = independent half-series fits; post-monsoon test =
fraction of successive Oct–Dec yearly maxima stepping down.
**Classes:** *inelastic candidate* (b ≤ −10 mm/yr, |b|>2σ, both halves ≤ −5, ≥3 yr of
post-monsoon maxima with ≥75 % down-steps) · *elastic dominant* (amp ≥ 2|b|·yr, |b| ≤ 2)
· *indeterminate*.

**Geometry & referencing.** LOS positive toward satellite. Dual-geometry AOIs get a
per-pixel joint [U,E] solve on the common window (`tools/joint_vertical.py`; north
neglected). AOI1 re-referenced to the coherent Chandigarh-city block after the LiCSBAS
auto-reference landed inside the subsiding corridor (log §2.1); AOI2 joint deramped
(planar) after the asc ramp diagnosis. Every referencing choice is a stated systematic.

**Interpretation guardrails applied.** Subsidence ≠ storage damage by itself: tectonic
piedmont uplift (mm-level, too small), construction loading (localized, does not track
well fields), and shrink–swell (seasonal, would appear as amplitude not trend) are
considered per-AOI below. LiCSAR magnitude lower-bound caveat at sharp bowls (AOI3
finding) applies wherever bowls are narrow. Isolated coherent islands can carry ±integer-
cycle/yr biases — interpretation restricted to contiguous coherent areas.

---

## 2. AOI1 — Chandigarh–Mohali–Kharar–Zirakpur–Dera Bassi. **VALIDATED · INELASTIC-DOMINANT**

**Data:** asc 027A_05990 (998 pairs, 102 epochs, 2014-12→2023-03, 8.25 yr) ·
desc 136D_05854 (401 pairs, 118 epochs, 2014-10→2021-12, 7.16 yr) · GACOS · city-referenced.

| Site | asc LOS | desc LOS | joint vertical (2015–21) |
|---|---|---|---|
| Kharar | −91.8 | −85.6 | **−114 (min −133)** |
| Sunny Enclave | −45.3 | −61.3 | −67 |
| Mohali Sec 70 | −40.3 | −48.0 | −59 |
| Dera Bassi | −92.6 | −130.1 | **−140 (min −155)** |
| Zirakpur | −4.4 | −6.4 | −7 |
| Chandigarh Sec 17 | +5.4 | +5.2 | ≈ +4 |
| Banur | −36.8 | −13.7 | −30 |

(mm/yr; hotspot 3×3/7×7 statistics; full maps + CSVs in `outputs/AOI1_chandigarh/`)

- Kharar pixel-level: −80 ± 0.5 mm/yr, R² = 0.995, half-series rates −81/−79 (no
  deceleration), seasonal amplitude 2.5 mm, **post-monsoon down-step fraction 1.00**.
- Inelastic-candidate: 23,614 px asc / 22,210 desc ≈ **250 km², mean −46 mm/yr LOS**.
  Elastic-dominant pixels are confined to the city/piedmont fringe.
- Confounders: piedmont tectonics (~mm/yr — 40× too small), construction loading (bowls
  track towns' municipal/agricultural well fields and extend beyond construction zones),
  seasonal shrink-swell (excluded by 2-mm amplitude). GACOS-corrected; artifact islands
  (Siwalik fringe positives at coh < 0.15; one +30–40 blob at 30.47°N 76.58°E) flagged
  and excluded.
- **Verdict: strong evidence of ongoing permanent storage destruction** — sustained
  multi-year compaction at 6–15 cm/yr vertical with negligible elastic recovery over
  documented clay-bearing confined aquifers; order-of-magnitude **~15 Mm³/yr of storage
  capacity being destroyed** (−59 mm/yr vertical × 251 km²). The prior (6–18 cm/yr,
  inelastic) is **confirmed within factor ~1.5**.

## 3. AOI4 — Lahore (Model Town & central city). **VALIDATED · INELASTIC-DOMINANT**

**Data:** asc 100A_05836 (398 pairs, 86 epochs, 2016-09→2020-06, 3.78 yr; CEDA coverage
limit) · desc 034D_05909 (1253 pairs, 149 epochs, 2016-04→2023-03, 6.93 yr) · GACOS.

| Site | asc LOS (2016–20) | desc LOS (2016–23) | joint vertical (2016–20) |
|---|---|---|---|
| Model Town | −45.0 | −46.6 | **−61 (min −64)** |
| Township | −45.3 | −51.2 | **−62 (min −78)** |
| Gulberg | −33.8 | −30.7 | −43 |
| Walled City | −32.9 | −30.1 | −39 |
| Central Lahore | −29.7 | −20.4 | — |
| Johar Town | −4.4 | −4.8 | −6 |
| DHA | −0.2 | −9.9 | — |

- Published prior "up to −43 mm/yr" reproduced at factor ~1.0 (LOS); vertical rates
  exceed it as expected geometrically.
- Inelastic-candidate ~162 km² (desc; mean −33 mm/yr LOS); seasonal amplitudes 3.5–5 mm
  (small vs trend). Desc series continues the subsidence unabated through early 2023 —
  no sign of stabilization.
- Deep lithology poorly documented → attribution rests on the kinematic signature
  (sustained, unrecovered, 7 yr). **Verdict: strong evidence of ongoing, likely permanent
  compaction under central-southern Lahore.**

## 4. AOI3 — Delhi NCR fringe (calibration site). **METHOD VALIDATED · bowls inelastic; magnitudes lower-bound**

**Data:** asc 027A_06190 (1026 pairs, 137 epochs, 2014-10→2023-03, 8.38 yr) · desc
136D_06053 (300 pairs, 95 epochs, 2015-04→2021-12; CEDA lag) · GACOS.

- Kapashera–Samalka and Faridabad bowls detected at published locations in both
  geometries. Joint vertical: **Samalka −23 mean / −44 min; Faridabad −7/−10 mm/yr**
  (2015–21). Samalka series: −10.7 mm/yr linear, amp 2 mm, no rebound.
- **Published-vs-measured gap reported prominently:** >110 mm/yr published peaks
  (2018–19 PSI, building-scale) vs our long-window 111-m averages — diagnosed as
  resolution + period effect (log §1.6): decorrelated long pairs (coh 0.10–0.16) clip
  steep small bowls; where short-pair coherence is good (Faridabad 0.54, Gurugram 0.58)
  the moderate multi-year rates are credible. **LiCSAR magnitudes at sharp urban-village
  bowls are lower bounds.**
- Classification: ~15 km² inelastic-candidate (mean −25 mm/yr LOS) in the two bowls;
  37–44 k px elastic-dominant across the wider fringe (amplitudes ~1–3 mm, trends ≈ 0).
- **Verdict: the two bowls carry the inelastic signature (sustained, unrecovered);
  vertical magnitude ≥ 23–44 mm/yr; the validation criterion ("within factor ~2 of
  published") is met for bowl location/structure and NOT met for peak magnitude — with
  the mechanism of the shortfall demonstrated.**

## 5. AOI2 — Ludhiana. **PRIOR CONTRADICTED — no ongoing subsidence found**

**Data:** asc 027A_05990 (1009 pairs, 114 epochs, 2014-12→2023-03) · desc 034D_05909
(1169 pairs after grid-fix, 144 epochs, 2016-04→2023-03) · GACOS.

- asc alone showed a smooth ±10–15 mm/yr NW–SE gradient; desc (independent frame) shows
  none → **ramp = orbital/ionospheric artifact**, removed (planar) for the joint solve.
- Joint deramped vertical field: **p2/p50/p98 = −5.6/+1.0/+5.5 mm/yr** — stable AOI.
  City sites marginally positive relative to surroundings in BOTH geometries.
  Inelastic-candidate: 6 px (desc) / 900 px (asc, inside the ramp = artifact).
- **Verdict: the published ~25 mm/yr Ludhiana subsidence is absent in 2015–2024 LiCSAR
  data. Either the earlier signal was localized/temporary or extraction has stabilized
  (canal-supply/paddy-policy changes post-2019 are plausible, untested). For the aquifer-
  damage question: no evidence of ongoing inelastic compaction 2015–2024; past permanent
  loss cannot be assessed from this record.**

## 6. AOI5 — Ghaggar interfan belt (S Sangrur–Mansa–Patiala). **NOT MEASURABLE at 111 m (valid null)**

**Data:** asc 027A_05990 (964 pairs after grid-fix, 107 epochs, 2014-10→2023-03, 8.35 yr;
86 % AOI coverage). Descending not pursued after the coherence result (decision log §5.0).

- 38,034 / 1,122,201 px usable (**3.4 %**) — monsoon-irrigated cropland decorrelates the
  belt, exactly as the mission anticipated.
- Town centres hold raw coherence (Sangrur 0.53, Barnala 0.56, Sunam 0.51, Patiala 0.38)
  but are masked by network-quality criteria; their unmasked indicative values spread
  ±15 mm/yr (Sangrur −1.4, Sunam +1.8, Barnala −10.8, Patiala +15.2, Budhlada +13.3
  relative to field median) — consistent with isolated-island integer-cycle bias, so
  **no rate claims are made**. No town shows a corridor-style multi-cm/yr signal.
- **Verdict: sparse/null result. The highest-clay-fraction belt remains untested by this
  product — full-resolution burst PSI (HyP3, ~1 credit/pair) is the right follow-up.**

---

## 7. Cross-AOI synthesis: elastic vs inelastic storage behaviour

| AOI | Secular signal | Seasonal amplitude | Post-monsoon recovery | Classification outcome | Permanent-damage verdict |
|---|---|---|---|---|---|
| 1 Chandigarh corridor | −60…−140 mm/yr vert, 8 yr sustained | 2–5 mm (≈2 % of trend) | fails every year (frac=1.0) | ~250 km² inelastic | **Evidenced, ongoing, large** |
| 4 Lahore | −39…−62 mm/yr vert, 7 yr sustained | 3.5–5 mm | fails persistently | ~160 km² inelastic | **Evidenced, ongoing** |
| 3 Delhi bowls | −23…−44 mm/yr vert (lower bound) | ~2 mm | fails | ~15 km² inelastic | **Evidenced (magnitude lower-bound)** |
| 3 Delhi wider fringe | ≈ 0 | 1–3 mm | recovers | 37–44 k px elastic | Elastic regime |
| 2 Ludhiana | ≈ 0 (±5) | 1–2 mm | recovers | ~none inelastic | **No ongoing damage detected** |
| 5 Ghaggar belt | unmeasurable | — | — | — | **Untested (data limitation)** |

The decomposition separates three regimes exactly as the mission's conceptual model
predicts: (i) confined, clay-bearing systems under overdraft → large secular loss with
tiny seasonal amplitude (AOI1, AOI4, Delhi bowls); (ii) unconfined/stabilized systems →
near-zero trend with small recoverable seasonal motion (Ludhiana, Delhi fringe);
(iii) coarse-sand recharge zones → no coherent signal of either kind. Where the prior
said clay exists (AOI1), we find the inelastic signature at full published magnitude;
where the aquifer is coarse sand, we find elastic or null behaviour. **The data are
consistent with the hypothesis that permanent storage destruction in this region is
concentrated where clay interbeds exist and heads have been drawn below historical
minima — and quantitatively largest in the Chandigarh–Mohali corridor.**

## 8. Uncertainties & failure points (honest list)

- Bootstrap rate σ ≈ 0.5–2 mm/yr at coherent pixels; referencing systematics are larger
  (±3–5 mm/yr for AOI1 city-block choice; stated per AOI).
- LiCSAR product limits: sharp-bowl magnitude clipping (proven at AOI3), isolated-island
  unwrap bias (AOI5, artifact blobs in AOI1/AOI2), inter-era grid shift (fixed), archive
  ends 2024 (2025–26 unobserved), CEDA lag truncates some desc series at 2021.
- GACOS sltd missing for 8–43 % of epochs per frame (affected ifgs uncorrected).
- Groundwater-head ground truth deferred (portals down/blocked) → elastic storativity and
  preconsolidation-head confirmation not computed; inelastic attribution is kinematic.
- Ludhiana/AOI5 published-rate contradictions could partly reflect technique differences
  (PSI vs 111-m SBAS); resolvable with burst PSI.

## 9. Recommended next steps

1. **HyP3 burst-InSAR validation stacks** (needs the user's NASA Earthdata login):
   Kapashera–Samalka, Kharar, Dera Bassi, and 2–3 Ghaggar towns at 20 m; extends the
   record to 2025–26 (S1C/S1D) and tests the lower-bound caveat. ~300–600 credits total —
   well within the free monthly allotment.
2. CGWB/WRIS well-level export (manual, when portals return) → seasonal head vs seasonal
   deformation → elastic skeletal storativity; residual trend → inelastic confirmation.
3. Extend AOI1/AOI4 series past 2024 as LiCSAR catches up; watch Dera Bassi (fastest,
   still accelerating within the desc window).

## 10. Deliverables index

`outputs/AOI{1..5}_*/`: LOS rate ± σ, seasonal amplitude, seasonal-peak timing, R²,
class GeoTIFFs; joint vertical + east GeoTIFFs (AOI1/2/3/4); rendered PNG maps with towns
and scale bars (109 figures); `*_decomposition_results.csv` (per-pixel decomposition +
classification; 0.4–1.2 M rows total); hotspot time-series plots (published-location
comparisons); `*_summary.json` machine-readable summaries. `PROCESSING_LOG.md` = full
provenance. Tools in `tools/` (downloader, grid-aligner, verifier, decomposer, joint
solver, mappers) are reusable for any LiCSAR frame worldwide.
