#!/usr/bin/env bash
# set -e
# sbatch run_lag_ridgeregression_vox_perm.sh
#name the job run_lag_ridgeregression_vox_perm and place it's output in a file named slurm-<jobid>.out
# allow 40 minutes to run (it should not take 40 minutes however)
# set partition to 'all' so it runs on any available node on the cluster
# sbatch --array=1-2 run_lag_ridgeregression_vox_perm.sh

#SBATCH -J 'lag_ridgeregression_vox_perm'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 2-23


module load matlab

# matlab -nosplash -nojvm -nodisplay -nodesktop -r "try;lag_ridgeregression_vox_perm($SLURM_ARRAY_TASK_ID);catch;exit;end;exit;"
matlab -nojvm -r "lag_ridgeregression_vox_perm($SLURM_ARRAY_JOB_ID)"



