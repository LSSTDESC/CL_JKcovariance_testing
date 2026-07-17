# Claude Progress

## 2026-07-17: Configurable low-memory mode for TreeCorr

- `cljkcov/compute_treecorr_cov.py` now reads `low_memory_flag` from the yaml config (`yd.get('low_memory_flag', True)`, so the default is True when the key is absent) and passes it as `low_mem` to both `ng.process` and `rg.process`.
- Added `low_memory_flag: True` to all config files:
  - `cljkcov/configs/smallsurvey.yml`
  - `cljkcov/configs/t19lite.yml`
  - `cljkcov/configs/m14_4_15_z0_508_0_574.yml`
  - `cljkcov/validations/randoms_comparison/test_treecorr_randoms.yml`
  - `cljkcov/validations/randoms_comparison/test_treecorr_norandoms.yml`
