_Authored by Kabelo Tsiane, 2026-10-20_

# Issues with the DESC data registry dataset `dp2_known_clusters` 0.1.0

This note records the known problems with the DP2 known-cluster catalog in the DESC data registry, how the cluster catalog for our DP2 cluster-lensing runs was derived from it, and what we changed. Numbers were measured on 2026-10-19.

## The dataset

| | |
|---|---|
| Registry id / name / version | 133, `dp2_known_clusters`, 0.1.0 (owner `dericj`), registered 2026-09-21 |
| File | `/global/cfs/cdirs/lsst/utilities/data-registry/lsst_desc_working/user/dericj/dp2_known_clusters/dp2_cluster_overlap_jul23.csv` (7.3 MB) |
| Upstream | https://github.com/LSSTDESC/DP2_known_clusters, "jul23" build (no commit recorded) |
| Content | 55,910 rows, one per catalogue entry; `dp2_tract >= 0` selects 7,330 rows inside a DP2 tract |
| Columns | `RA, DEC, redshift, redshift_type, redshift_unc, bcg_specz, richness, richness_unc, M200m, M500c, catalog_source, catalog_id, dp2_tract` |

Rows per `catalog_source` (all / in a DP2 tract):

| catalogue | all | in DP2 |
|---|---:|---:|
| DES Y3 RedMaPPer | 21,092 | 2,759 |
| eRASS 1 | 12,247 | 517 |
| ACT DR6 | 10,040 | 1,661 |
| HSC Wide | 7,939 | 1,868 |
| MCXC-II | 2,221 | 233 |
| Planck | 1,653 | 171 |
| SPT ECS | 470 | 15 |
| HSC Deep | 248 | 106 |

## Issues in the registry dataset

1. **`catalog_id` is not a unique identifier.**
   - `(catalog_source, catalog_id)` is unique, but `catalog_id` alone is not: 9,696 repeated values overall, 473 among the DP2 rows. For example, DES and HSC both have small integer ids.
   - For MCXC-II (0–2220), Planck (1–1653), HSC Wide (1–7939) and HSC Deep (1–248) the ids are row numbers, not the catalogues' own ids, so rows cannot be traced back to the source catalogues.
   - The column mixes integers (DES, HSC, ...) and names (ACT, SPT, eRASS), so it is read as text.
   - Downstream code that matches clusters by `catalog_id` alone silently mixes clusters (see the derived file below).
2. **Cross-catalogue duplicates are not flagged.** The description says the file "contains cross-catalogue duplicates", but there is no flag or group column.
   - Among the DP2 rows there are 1,099 pairs from different catalogues within 1′ and |Δz| < 0.05. They merge into 6,373 groups from 7,330 rows (about 13% redundancy).
   - Largest overlaps: DES–HSC Wide 284, ACT–DES 219, DES–eRASS 121, ACT–HSC Wide 97, ACT–eRASS 96, ACT–MCXC 61, HSC Deep–HSC Wide 56, ACT–Planck 28.
3. **Richness definitions are not documented.**
   - DES richness starts at 20.0 (redMaPPer λ ≥ 20).
   - HSC Wide starts at exactly 15.0 and HSC Deep at 15.02. This is probably CAMIRA N_mem, which is on a different scale from λ, but the file does not say.
   - Richness is -1 for ACT, eRASS, MCXC and Planck. SPT runs from 0 to 204.
   - Combining richness across catalogues is therefore not meaningful without this information.
4. **Missing values use -1**, including in the text column `redshift_type`:
   - `redshift` = -1 for 559 Planck rows (57 in DP2);
   - `M500c` = -1 for all DES and HSC rows and 1,805 eRASS rows;
   - `M200m` = -1 everywhere except ACT and DES;
   - `redshift_unc` and `richness_unc` = -1 for all HSC Wide/Deep, MCXC and Planck rows and 155 eRASS rows;
   - `bcg_specz` is only filled for DES (2,411), HSC Wide/Deep (3,032 + 135) and SPT (470).
5. **SPT uses 0 as a second missing-value marker:** 172 rows with richness = 0 (5 in DP2) and 22 with redshift = 0 (none in DP2).
6. **Mixed mass units.** eRASS `M500c` is in units of 1e13 M☉ (median 19.5), while every other catalogue uses M☉. This is documented only in the eRASS pointer dataset (132), not in 133.
7. **`redshift_type` has no common vocabulary:**
   - DES, HSC, Planck and SPT use "-1";
   - ACT uses `phot` / `spec`;
   - MCXC uses S / SP / P / E / U;
   - eRASS uses `photo_z` / `spec_z_boot` / `cg_spec_z` / `lit_z`.
8. **Provenance is incomplete.**
   - The registry execution (130) is a "Fabricated execution for dataset" with no site, configuration or start time.
   - The only declared input is the eRASS1 pointer 0.1.0 (dataset 132), which has since been superseded by 1.0.0 (dataset 134).
   - The other six catalogues, their versions, and the code that assigns `dp2_tract` are not recorded.
   - The description says "8 external catalogues" but lists 7 names, because HSC Wide and HSC Deep are counted separately.
9. **"In a DP2 tract" is not "in the shear footprint".** `dp2_tract >= 0` only means the centre lies in a DP2 tract, and DP2 coverage of those tracts is patchy.
   - Only about 16% of the cluster centres in our first derived catalog (633 of 4,032) fall inside the Rubin wide-field shear mask (registry 151).
   - Of the 2,759 in-DP2 redMaPPer clusters, 708 have Rubin metaDetect ΔΣ background sources.

## Our derived catalog `external_0723_x_dp2` (replaced)

Built on 2026-09-21 by `~/make_dp2_cat.py` (also cells 3–8 of `~/lsdb.ipynb`) and converted to HDF5 by what is now `~/txp/dp2/make_cluster_catalog_external_0723.py`:
- read a copy of the registry CSV, byte-identical to the registry file;
- matched each DP2 object (`object_collection` HATS table) to its nearest cluster within 1″;
- kept rows with redshift and richness ≠ -1;
- wrote `/pscratch/sd/k/kabelo/dp2/external_0723_x_dp2.{csv,hdf5}` (4,032 rows).

Problems found:
- **Mixed richness:** DES redMaPPer (2,385), HSC Wide (1,566) and HSC Deep (81) were combined in one `richness` column, and the same cluster could appear twice (65 DES–HSC or HSC Wide–Deep pairs within 1′, |Δz| < 0.05).
- **Ids not unique:** `cluster_id` = `catalog_id`, giving only 3,814 unique ids in 4,032 rows.
  - 194 ids were shared between DES and HSC (397 rows).
  - 17 clusters appeared twice, matched to two DP2 objects.
  - TXPipe's `CLClusterEnsembleProfiles` looked clusters up by id and took the first match, so a cluster sharing an id was given **another cluster's sources**. The stage now stops with an error if ids are not unique.
- **Clusters lost in the match:** the 1″ match dropped about 14% of in-DP2 clusters (DES 389 of 2,759, HSC Wide 304, HSC Deep 25, SPT all 15). The loss is roughly independent of richness and redshift.
- **Positions:** `ra, dec` were the matched DP2 object's position (0.12″ median offset from the catalogue centre), and the DP2 objectId was not kept.
- **Missing values passed on:** 1,647 HSC rows had `redshift_err` = `richness_err` = -1, which were passed on as real values.

## What we use now

`~/txp/dp2/make_cluster_catalog.py` builds `/pscratch/sd/k/kabelo/dp2/dp2_known_clusters_0.1.0_redmapper.hdf5` directly from registry 133:
- **Selection:** `catalog_source == "DES Y3 RedMaPPer"` and `dp2_tract >= 0`, giving 2,759 clusters with z 0.10–0.91 and λ 20–178.
- **Positions:** the catalogue's original `RA, DEC`.
- **Ids:** `cluster_id` = `catalog_id` as int64, checked to be unique.
- **Checks:** the build stops if there are -1 values in redshift, richness or their errors.
- **Provenance:** stored as HDF5 attributes (registry id, name, version, path, selection, build date).

## Suggestions for the dataset maintainers

- **Ids:** add a globally unique id column, and keep each source catalogue's own id.
- **Duplicates:** add a duplicate-group column, or a preferred-entry flag.
- **Richness:** document the richness definition and catalogue version for each `catalog_source`.
- **Missing values:** use empty values (NaN) or one documented marker per column, instead of -1 and 0.
- **Units and vocabulary:** give all masses in the same unit, and one vocabulary for `redshift_type`.
- **Provenance:** record the inputs (all catalogues and versions), the upstream commit, and the code that assigns `dp2_tract`.
- **Footprint:** optionally add a flag for clusters inside the DP2 shear footprint.
