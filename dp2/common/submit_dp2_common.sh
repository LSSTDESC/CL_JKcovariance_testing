#!/bin/bash
#SBATCH -A m1727
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00
#SBATCH --job-name=dp2_common
#SBATCH --output=/pscratch/sd/k/kabelo/dp2/outputs/slurm-%x-%j.out
#SBATCH --mail-user=kabelo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#
# Usage, from ~/txp/dp2:   sbatch -J dp2_common_<catalog> common/submit_dp2_common.sh <catalog>
#   <catalog> = rubin | sheldon   (runs common/pipeline_<catalog>.yml)
#
# Time limit: the Rubin run from scratch took ~35 min (Sep 2026), dominated by
# BPZliteEstimator (~17 min). The Sheldon catalog has ~2.8x as many rows, and
# TXCutShearCatalog and BPZliteEstimator scale with rows, so allow ~2.5 h.
# resume: true skips finished stages on reruns.

set -e
CATALOG=${1:?usage: sbatch common/submit_dp2_common.sh <rubin|sheldon>}
cd /global/u2/k/kabelo/CL_JKcovariance_testing/dp2/common
PIPELINE=pipeline_${CATALOG}.yml
[ -f "$PIPELINE" ] || { echo "No $PIPELINE in $(pwd)"; exit 1; }

# Uses the TXPipe conda env built by TXPipe/bin/perlmutter-install.sh. Its
# h5py/hdf5 use conda-forge's "external" mpich, and its activate.d hook
# (activate-txpipe.sh) loads cray-mpich and puts $MPICH_DIR/lib-abi-mpich on
# LD_LIBRARY_PATH, so h5py and mpi4py both link to the same Cray MPICH.
# Activating also puts the env's bin first on PATH, so the
# "python3 -m txpipe ..." subprocesses that ceci launches use the right
# interpreter.
TXPIPE_ENV=/global/cfs/cdirs/m1727/kabelo/TXPipe/conda
module load python
set +e  # conda activate hooks are not always "set -e" safe
conda activate $TXPIPE_ENV
set -e

CECI=$TXPIPE_ENV/bin/ceci

mkdir -p /pscratch/sd/k/kabelo/dp2/outputs/${CATALOG}/logs

echo "=== Starting DP2 shared-stage pipeline for ${CATALOG}: $(date) ==="
echo "TXPipe branch: $(git -C /global/cfs/cdirs/m1727/kabelo/TXPipe rev-parse --abbrev-ref HEAD) @ $(git -C /global/cfs/cdirs/m1727/kabelo/TXPipe rev-parse --short HEAD)"
$CECI $PIPELINE

echo "=== Done: $(date) ==="
