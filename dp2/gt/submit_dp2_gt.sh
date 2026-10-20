#!/bin/bash
#SBATCH -A m1727
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --job-name=dp2_gt
#SBATCH --output=/pscratch/sd/k/kabelo/dp2/outputs/slurm-%x-%j.out
#SBATCH --mail-user=kabelo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#
# Usage, from ~/txp/dp2:   sbatch -J dp2_gt_<catalog> gt/submit_dp2_gt.sh <catalog>
#   <catalog> = rubin | sheldon   (runs gt/pipeline_<catalog>.yml; run common/pipeline_<catalog>.yml first)
#
# Time limit: only the two cluster stages run here. Jobs took 5-6 min (Sep/Oct
# 2026) but up to 13 min on a slow filesystem day, so allow 30 min.
# resume: true skips finished stages on reruns.

set -e
CATALOG=${1:?usage: sbatch gt/submit_dp2_gt.sh <rubin|sheldon>}
cd /global/u2/k/kabelo/CL_JKcovariance_testing/dp2/gt
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

mkdir -p /pscratch/sd/k/kabelo/dp2/outputs/${CATALOG}/gt/logs

echo "=== Starting DP2 gt pipeline for ${CATALOG}: $(date) ==="
echo "TXPipe branch: $(git -C /global/cfs/cdirs/m1727/kabelo/TXPipe rev-parse --abbrev-ref HEAD) @ $(git -C /global/cfs/cdirs/m1727/kabelo/TXPipe rev-parse --short HEAD)"
$CECI $PIPELINE

echo "=== Done: $(date) ==="
