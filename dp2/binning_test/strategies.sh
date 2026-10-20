# Cluster binnings compared in binning_test (sourced by submit_binning_test.sh).
# Catalog: dp2_known_clusters 0.1.0, DES Y3 redMaPPer in DP2 tracts (lambda >= 20).
# Weakest bin, as clusters with sources : occupied k-means regions (30 shared regions,
# seed 11, from the rubin delta_sigma cluster_shear_catalogs of 2026-10-19):
#   A 37:21, B 37:21, C 91:28, D 36:19 (the lambda 45-200 bin).
declare -A ZEDGE=(
    [A]="[0.2, 0.35, 0.5, 0.8]"   # baseline, same as common/config.yml
    [B]="[0.2, 0.35, 0.5, 0.8]"   # richness-edge sensitivity
    [C]="[0.2, 0.45, 0.8]"        # 2 vs 3 redshift bins
    [D]="[0.2, 0.45, 0.8]"        # adds a high-richness bin
)
declare -A RICHEDGE=(
    [A]="[20., 24., 32., 200.]"        # 9 bins
    [B]="[20., 25., 35., 200.]"        # 9 bins
    [C]="[20., 24., 32., 200.]"        # 6 bins
    [D]="[20., 24., 32., 45., 200.]"   # 8 bins
)
