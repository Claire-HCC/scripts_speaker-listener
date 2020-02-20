#!/bin/bash
#SBATCH -J 'matlab_job'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 3:0:00
#SBATCH -c 10
#SBATCH --array=1-1000
#SBATCH --mem-per-cpu=2G


module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

perm=$SLURM_ARRAY_TASK_ID
declare -i perm


# matlab -nodisplay -r "roi_tr_pattern_regression_SL_lagSelection($perm)"
# matlab -nodisplay -r "roi_tr_pattern_regression_LL_lagSelection($perm)"
# matlab -nodisplay -r "roi_tr_pattern_regression_SL_permSubj($perm)"
# matlab -nodisplay -r "roi_tr_pattern_regression_LL_permSubj($perm)"

matlab -nodisplay -r "wholeBrain_tr_temporal_lagcorr_SLg_permPhase($perm)"


# matlab -nodisplay -r "roi_tr_bined_pattern_granger_SL_permSubj($s)"
# matlab -nodisplay -r "roi_tr_bined_pattern_regression_LL_permSubj($s)"









