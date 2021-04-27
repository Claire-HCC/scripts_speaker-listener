#!/bin/bash
#SBATCH -J 'subj2'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 3:00:00
#SBATCH -c 5
#SBATCH --array=1-18
#SBATCH --mem-per-cpu=500M


module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

perm=$SLURM_ARRAY_TASK_ID
declare -i perm

matlab -nodisplay -r "roi_bined_pattern_regression_SL_permSubj($perm)"
# matlab -nodisplay -r "roi_bined_pattern_regression_SL_permSubj($perm)"

# sqeueu -u huichuan
# scancel <jobid>









