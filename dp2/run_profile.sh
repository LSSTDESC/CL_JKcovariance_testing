#!/bin/bash
LOGDIR=/pscratch/sd/k/kabelo/dp2/outputs

for catalog in ${@:-rubin sheldon}; do
    for profile in delta_sigma gt; do
        log=$LOGDIR/${profile}_${catalog}.log
        echo "##### $profile $catalog -> $log"
        bash $profile/submit_dp2_${profile}.sh $catalog 2>&1 | tee $log
    done
done
echo "##### all done: $(date)"
