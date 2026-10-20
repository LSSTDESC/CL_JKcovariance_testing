#!/bin/bash
#SBATCH -A m1727
#SBATCH -C cpu
#SBATCH -q shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=00:45:00
#SBATCH --job-name=takahashi_jk_compare_100deg2
#SBATCH --output=/pscratch/sd/k/kabelo/txp_cluster_compare/logs_t17/slurm-%j.out
#SBATCH --mail-user=kabelo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -e
cd /global/u2/k/kabelo/CL_JKcovariance_testing

module load cray-mpich-abi

# ceci launches each stage as a subprocess via "python3 -m txpipe ...", which
# otherwise resolves to the system python3.6 (no ceci/txpipe installed) rather
# than the txpipe_clp venv -- put the venv's bin first so those subprocesses
# pick up the right interpreter too.
export PATH=/global/cfs/cdirs/m1727/kabelo/TXPipe/txpipe_clp/bin:$PATH

CECI=/global/cfs/cdirs/m1727/kabelo/TXPipe/txpipe_clp/bin/ceci

echo "=== Starting Takahashi healpix vs kmeans jackknife g_t comparison pipeline: $(date) ==="
$CECI takahashi/jackknife/takahashi_jk_compare_100deg2.yml

echo "=== Done: $(date) ==="
