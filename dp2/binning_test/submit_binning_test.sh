#!/bin/bash
#SBATCH -A m1727
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --job-name=dp2_binning_test
#SBATCH --output=/pscratch/sd/k/kabelo/dp2/outputs/slurm-%x-%j.out
#SBATCH --mail-user=kabelo@umich.edu
#SBATCH --mail-type=FAIL
#
# Usage, from ~/txp/dp2:
#   sbatch -J dp2_bin<strategy>_<profile>_<catalog> binning_test/submit_binning_test.sh <strategy> <catalog> <profile>
#   <strategy> = A | B | C | D (strategies.sh), <catalog> = rubin | sheldon, <profile> = delta_sigma | gt
# or submit all of them (16 with A-D) with binning_test/submit_all.sh. Needs outputs/<catalog>/<profile>/cluster_shear_catalogs.hdf5.
#
# Time limit: the binning stage takes seconds and CLClusterEnsembleProfiles took
# ~1-4 min (2026-10-19/20 runs), so 30 min is plenty.

set -e
STRATEGY=${1:?usage: submit_binning_test.sh <A|B|C|D> <rubin|sheldon> <delta_sigma|gt>}
CATALOG=${2:?usage: submit_binning_test.sh <A|B|C|D> <rubin|sheldon> <delta_sigma|gt>}
PROFILE=${3:?usage: submit_binning_test.sh <A|B|C|D> <rubin|sheldon> <delta_sigma|gt>}
cd /global/u2/k/kabelo/CL_JKcovariance_testing/dp2/binning_test
source strategies.sh
[ -n "${RICHEDGE[$STRATEGY]}" ] || { echo "Unknown strategy $STRATEGY"; exit 1; }
[ -f ../$PROFILE/config.yml ] || { echo "No ../$PROFILE/config.yml"; exit 1; }

OUT=/pscratch/sd/k/kabelo/dp2/outputs/${CATALOG}/binning_${STRATEGY}/${PROFILE}
mkdir -p $OUT/logs

# Profile settings come from ../<profile>/config.yml unchanged; only the
# binning stage section is added. A copy is kept next to the outputs.
CONFIG=$OUT/config.yml
cat ../$PROFILE/config.yml > $CONFIG
cat >> $CONFIG <<EOF

# Added by binning_test/submit_binning_test.sh, strategy $STRATEGY
CLClusterBinningRedshiftRichness:
    zedge: ${ZEDGE[$STRATEGY]}
    richedge: ${RICHEDGE[$STRATEGY]}
EOF

# Same TXPipe conda env setup as ../delta_sigma/submit_dp2_delta_sigma.sh
TXPIPE_ENV=/global/cfs/cdirs/m1727/kabelo/TXPipe/conda
module load python
set +e  # conda activate hooks are not always "set -e" safe
conda activate $TXPIPE_ENV
set -e

echo "=== Starting DP2 binning test $STRATEGY for $CATALOG $PROFILE: $(date) ==="
echo "richedge: ${RICHEDGE[$STRATEGY]}  zedge: ${ZEDGE[$STRATEGY]}"
echo "TXPipe branch: $(git -C /global/cfs/cdirs/m1727/kabelo/TXPipe rev-parse --abbrev-ref HEAD) @ $(git -C /global/cfs/cdirs/m1727/kabelo/TXPipe rev-parse --short HEAD)"
$TXPIPE_ENV/bin/ceci pipeline.yml -t catalog=$CATALOG profile=$PROFILE strategy=$STRATEGY config=$CONFIG ${DRY_RUN:+--dry-run}

echo "=== Done: $(date) ==="
