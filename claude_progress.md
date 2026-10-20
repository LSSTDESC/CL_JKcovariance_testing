# Claude progress log

## 2026-10-20: Imported DP2 / Takahashi jackknife covariance testing from ~/txp

Copied (originals left in place) from `~/txp`, which is not a git repo:
- `dp2/`: from `~/txp/dp2`. This is the DP2 TXPipe pipelines comparing kmeans and healpix jackknife covariance for ΔΣ and g_t.
  - `common/`, `delta_sigma/`, `gt/`: current `config.yml`, `pipeline_{rubin,sheldon}.yml`, `submit_*.sh` and flowchart PNGs.
  - `binning_test/`: JK-region occupancy per (z, richness) bin, and binning strategies A–D.
  - top level: `DP2_CL_pipeline.ipynb`, `viz.ipynb`, `run_profile.sh`, `submit_dp2_chain.sh`, `make_cluster_catalog*.py`, `REGISTRY_dp2_known_clusters_issues.md`.
- `takahashi/`:
  - `jackknife/`: the healpix vs kmeans g_t comparison pipeline, with `submit_takahashi_jk_compare_100deg2.sh`.
  - `make_t17_cutout.py` and `validate_t17_shear_convention.py`.
- `cosmology/cosmology_t17_fid.yml`: shared by all the pipelines above.

Not copied:
- `.pre_*` and `.bak` backups, `.ipynb_checkpoints/`
- `rail-bpz-inputs/` and `bpz_riz.columns`
- the Confluence HTML
- `dp2_shear_catalog_comparison.ipynb`
- `~/txp/{dc2,deprecated,source_select_outputs}`
- all `/pscratch` data

Path rewrites, the only content changes:
- `~/txp/dp2` → `CL_JKcovariance_testing/dp2`. This covers the `cd` lines in the submit scripts, `my_dp2_dir` in the notebook, `viz.ipynb`, and the `exec(open(...jk_regions.py))` in `jk_fixed.py` and `jk_fullcat.py`.
- `~/txp/cosmology_t17_fid.yml` → `CL_JKcovariance_testing/cosmology/cosmology_t17_fid.yml`.
- Takahashi:
  - The submit script now does `cd`s to the repo root.
  - The pipeline uses `./cosmology/cosmology_t17_fid.yml`.

`/pscratch` paths, the TXPipe checkout (`/global/cfs/cdirs/m1727/kabelo/TXPipe`, branch `cl_kmeans_jackknife`) and its env are unchanged.

Checks run:
- `diff` against the originals shows only the rewrites.
- No `kabelo/txp` references remain.
- `bash -n`, `py_compile`, and notebook JSON / YAML parsing all pass.
- ceci `--dry-run` was not completed: h5py in the txpipe_clp venv can't find `libmpi.so.12` in a login shell.

Known issues, carried over unfixed:
- `DP2_CL_pipeline.ipynb` §5.3 does `git show cl_kmeans_jackknife:.../jk_kmeans.py`. That file was deleted in TXPipe ef57b26e; ref `b6545049` still has it.
- The Takahashi pipeline references the TXPipe module `plot_jackknife_comparison` (`CLClusterJackknifeComparisonPlots`), which no longer exists on the branch.
- `dp2/binning_test/jk_regions.py` reads the old `external_0723_x_dp2.hdf5` catalog, and its bins differ from `strategies.sh`.
