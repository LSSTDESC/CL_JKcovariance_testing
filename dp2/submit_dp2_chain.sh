#!/bin/bash
# Submit the full DP2 cluster-lensing chain for one shear catalog:
#   common  ->  delta_sigma and gt (in parallel, each only if common succeeds)
#
# Usage, from ~/txp/dp2:   ./submit_dp2_chain.sh <catalog> [--skip-common]
#   <catalog> = rubin | sheldon
#   --skip-common: common has already run; submit only the two profile pipelines
#
# SLURM logs: /pscratch/sd/k/kabelo/dp2/outputs/slurm-<jobname>-<jobid>.out

set -e
CATALOG=${1:?usage: ./submit_dp2_chain.sh <rubin|sheldon> [--skip-common]}
cd "$(dirname "$0")"
for d in common delta_sigma gt; do
    [ -f "$d/pipeline_${CATALOG}.yml" ] || { echo "No $d/pipeline_${CATALOG}.yml"; exit 1; }
done

if [ "$2" == "--skip-common" ]; then
    DEP=""
else
    jid=$(sbatch --parsable -J dp2_common_${CATALOG} common/submit_dp2_common.sh ${CATALOG})
    echo "common:      job $jid"
    DEP="--dependency=afterok:$jid"
fi
echo "delta_sigma: job $(sbatch --parsable $DEP -J dp2_delta_sigma_${CATALOG} delta_sigma/submit_dp2_delta_sigma.sh ${CATALOG})"
echo "gt:          job $(sbatch --parsable $DEP -J dp2_gt_${CATALOG} gt/submit_dp2_gt.sh ${CATALOG})"
