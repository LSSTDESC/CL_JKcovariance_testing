#!/bin/bash
# Submit every binning-test combination: strategies x catalogs x profiles.
#
# Usage, from ~/txp/dp2:   binning_test/submit_all.sh [strategy ...]   (default: A B C D)
#
# SLURM logs: /pscratch/sd/k/kabelo/dp2/outputs/slurm-dp2_bin<strategy>_<profile>_<catalog>-<jobid>.out

set -e
cd "$(dirname "$0")/.."
for s in ${@:-A B C D}; do
    for c in rubin sheldon; do
        for p in delta_sigma gt; do
            echo "$s $c $p: job $(sbatch --parsable -J dp2_bin${s}_${p}_${c} binning_test/submit_binning_test.sh $s $c $p)"
        done
    done
done
