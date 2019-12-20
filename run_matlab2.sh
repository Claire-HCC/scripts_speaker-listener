#!/bin/bash
#SBATCH -J 'matlab_job'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 3:00:00
#SBATCH -c 5
#SBATCH --array=1-2
#SBATCH --mem-per-cpu=500M


module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# perm=$SLURM_ARRAY_TASK_ID
# declare -i perm

matlab -nodisplay -r "roi_tr_pattern_regression_SL"
# matlab -nodisplay -r "roi_tr_pattern_regression_LL"









