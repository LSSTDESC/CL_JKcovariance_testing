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

## 2026-10-20 / 2026-11-03: DP2 Sheldon profiles, 0.2–1.5 Mpc, k-means-30 jackknife (branch `dp2_1halo`)

Profile settings (`dp2/{delta_sigma,gt}/config.yml`, shared by Rubin and Sheldon):
- ΔΣ: 5 log bins over 0.2–1.5 Mpc (edges 0.2, 0.299, 0.448, 0.670, 1.002, 1.5). Sources selected within `max_radius` 1.8 Mpc (was 5).
- g_t: 5 log bins over 0.6–4.4 arcmin, which is 0.2–1.5 Mpc at z≈0.43. Sources selected within `max_angle` 5 arcmin (was 22).
- Jackknife unchanged: `jackknife_regions: kmeans`, 30 regions, seed 11.

Output paths: the profile pipelines now write to `outputs/<cat>/{delta_sigma,gt}_1p5mpc/`. The old 6-bin outputs are still in `outputs/<cat>/{delta_sigma,gt}/`.

Run: Sheldon only, in an interactive job on 2026-10-20 (`./run_profile.sh sheldon`), TXPipe `cl_kmeans_jackknife` @ ef57b26e.
- ΔΣ: about 4 min in total.
- g_t: about 2.5 min.
- Regions: 30 k-means regions from about 610 clusters with sources. Each bin uses 18–28 of them.

`dp2/DP2_CL_pipeline.ipynb` was rewritten for this run and executed with its outputs saved. Sections:
1. settings and how to run
2. profiles and correlation matrices
3. k-means regions on the sky
4. k-means vs healpix nside 16, error against radius

Findings:
- σ_kmeans/σ_healpix: median 0.97, consistent with jackknife noise.
- Cross component inconsistent with zero in 3 of 9 bins, the same bins for ΔΣ and g_t. To follow up.

Branches: this work was first done on `clmm_testing`, then moved to `dp2_1halo`, which starts from `main`. The unpushed cherry-pick of the dp2 import (1f0dc5a) was removed from local `clmm_testing`, which now matches origin again (d360cb5).
