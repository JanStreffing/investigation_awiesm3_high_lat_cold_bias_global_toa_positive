#!/bin/bash -l
#SBATCH --job-name=submit_PI200_branches
#SBATCH --account=ab0246
#SBATCH --partition=shared
#SBATCH --ntasks=4
#SBATCH --time=02:00:00
#SBATCH --output=/work/bb1469/a270092/runtime/logs_branch_submit/submit_PI200_branches_%j.log
# Submit the four 1560 branches of PI200 once its 1550-1559 leg has written the 1560
# restarts.  Launch with --dependency=afterok:<PI200 leg job>.  The check-mode (-c)
# directories are moved aside first: an existing experiment directory makes
# esm_runscripts prompt, and it staged no restarts because they did not exist yet.
set -u
# sbatch exports this job's environment to the jobs it submits. A SLURM_MEM_PER_NODE left
# over from here collides with the compute partition's SLURM_MEM_PER_CPU and srun refuses to
# start ("mutually exclusive"), after which esm_tools resubmits through every leg.
for v in $(env | grep -o '^SLURM[A-Z_]*'); do unset $v; done
B=/work/bb1469/a270092/runtime/awiesm3-v3.4
R=$B/PI200/restart
RS=/home/a/a270092/esm_tools/runscripts/awiesm3/develop
STAMP=$(date +%Y%m%d_%H%M)

n_oce=$(ls $R/fesom/fesom.1560-01-01.oce.restart 2>/dev/null | wc -l)
n_ice=$(ls $R/fesom/fesom.1560-01-01.ice.restart 2>/dev/null | wc -l)
n_lpj=$(ls $R/lpj_guess/lpjg_state_1560 2>/dev/null | wc -l)
ref_oce=$(ls $R/fesom/fesom.1550-01-01.oce.restart | wc -l)
ref_ice=$(ls $R/fesom/fesom.1550-01-01.ice.restart | wc -l)
ref_lpj=$(ls $R/lpj_guess/lpjg_state_1550 | wc -l)
echo "1560 restarts: oce $n_oce/$ref_oce ice $n_ice/$ref_ice lpjg $n_lpj/$ref_lpj"
if [ "$n_oce" -ne "$ref_oce" ] || [ "$n_ice" -ne "$ref_ice" ] || [ "$n_lpj" -ne "$ref_lpj" ] || [ "$n_oce" -eq 0 ]; then
    echo "1560 restarts incomplete, nothing submitted"; exit 1
fi

cd $RS
for a in spp mle h0 all3; do
    E=$B/PI200_$a
    if [ -d $E ]; then
        if ls $E/outdata/*/* >/dev/null 2>&1; then echo "$E holds output, not touching it"; continue; fi
        mv $E $E.checkmode_$STAMP
    fi
    echo "== submitting PI200_$a"
    esm_runscripts awiesm3-develop-levante-TCO95L91-CORE3_PI200_$a.yaml -e PI200_$a < /dev/null \
        > /work/bb1469/a270092/runtime/logs_branch_submit/PI200_${a}_$STAMP.log 2>&1
    echo "   exit $?"
    grep -i "Submitted batch job" /work/bb1469/a270092/runtime/logs_branch_submit/PI200_${a}_$STAMP.log | tail -1
    W=$E/run_15600101-15691231/work
    echo "   staged restarts: oce $(ls $W/fesom.1560-01-01.oce.restart 2>/dev/null | wc -l)" \
         "ice $(ls $W/fesom.1560-01-01.ice.restart 2>/dev/null | wc -l)" \
         "lpjg $(ls -d $W/lpjg_state_* 2>/dev/null | head -1)"
done
squeue -u $USER -o "%.10i %.14j %.8T %.20E"
