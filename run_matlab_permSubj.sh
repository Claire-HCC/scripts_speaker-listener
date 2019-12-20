#!/bin/bash
#SBATCH -J 'matlab_job'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 03:00:00
#SBATCH -c 15
#SBATCH --array=1-48
#SBATCH --mem-per-cpu=100M


module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

perm=$SLURM_ARRAY_TASK_ID
declare -i perm


matlab -nodisplay -r "roi_tr_pattern_regression_SL_permSubj($perm)"
matlab -nodisplay -r "roi_tr_pattern_regression_LL_permSubj($perm)"

# matlab -nodisplay -r "roi_tr_pattern_regression_LL_permPhase($perm)"


# matlab -nodisplay -r "roi_tr_bined_pattern_granger_SL_permSubj($s)"
# matlab -nodisplay -r "roi_tr_bined_pattern_regression_LL_permSubj($s)"









