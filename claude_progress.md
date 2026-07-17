# Claude Progress

## 2026-07-17: Area-based survey footprint selection (branch small_area_testing)

- Replaced the `nsky` longitude-wedge selection in `cljkcov/compute_treecorr_cov.py` with an area-based survey footprint: the config now specifies `survey_area` (deg^2), and a new function `select_survey_pixels` picks an ~circular list of healpixels (nside 64, pixel ~0.84 deg^2) around theta, phi = pi/2, pi (center of the coordinate ranges) using the spherical-cap radius + `hp.query_disc`. Selected area matches the request to within a few deg^2.
- Halos, shear map pixels, and randoms are all selected by membership of their nside-64 healpixel in the footprint; centering at phi = pi keeps the footprint away from the phi = 0/2pi boundary, so no phi wrapping is needed.
- Randoms (10x the cluster count) are drawn uniformly within the footprint by `sample_survey_randoms`: a random coarse footprint pixel + a random NESTED child healpixel at nside 2^17 (~1.6 arcsec centers). Benchmarked ~3x faster than the earlier bounding-box rejection sampling (`tests/benchmark_randoms.py`); also fixes the old non-uniform-in-theta draw.
- Tests: `tests/test_survey_selection.py` (footprint) and `tests/test_survey_randoms.py` (randoms containment + uniformity).
- Output filenames changed: `datavector_<sim>_area<A>deg2.txt` / `cov_<sim>_area<A>deg2.txt` (was `..._patchN<n>_<j>.txt`).
- Only `configs/smallsurvey.yml` converted (`survey_area: [1300, 650, 320]`, placeholder values; `njk_patch` keyed by area). The other 4 configs still use `nsky` and need conversion before running with the new code.
- `Covariance_plotting_Comparison.ipynb` converted to the survey_area scheme: reads `yaml_data['survey_area']`, loads `*_area<A>deg2.txt` files (one per sim per area, no more per-patch jj loop), panels annotated with the survey area. `validations/randoms_comparison/compare_treecorr_randoms.ipynb` still uses the old `patchN*`/`nsky` scheme.
- Environment note: use DESC python (`source /global/common/software/lsst/common/miniconda/setup_current_python.sh`); system python lacks healpy/treecorr.

## 2026-07-17: Configurable low-memory mode for TreeCorr

- `cljkcov/compute_treecorr_cov.py` now reads `low_memory_flag` from the yaml config (`yd.get('low_memory_flag', True)`, so the default is True when the key is absent) and passes it as `low_mem` to both `ng.process` and `rg.process`.
- Added `low_memory_flag: True` to all config files:
  - `cljkcov/configs/smallsurvey.yml`
  - `cljkcov/configs/t19lite.yml`
  - `cljkcov/configs/m14_4_15_z0_508_0_574.yml`
  - `cljkcov/validations/randoms_comparison/test_treecorr_randoms.yml`
  - `cljkcov/validations/randoms_comparison/test_treecorr_norandoms.yml`
