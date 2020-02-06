#!/bin/bash
#SBATCH -J 'matlab_exp'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 3:00:00
#SBATCH -c 10
#SBATCH --array=1-4
#SBATCH --mem-per-cpu=2G


module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

ei=$SLURM_ARRAY_TASK_ID
declare -i ei


# matlab -nodisplay -r "roi_tr_pattern_regression_SL_permSubj($perm)"
# matlab -nodisplay -r "roi_tr_pattern_regression_LL_permSubj($perm)"


# matlab -nodisplay -r "roi_tr_pattern_regression_LL($ei)"
# matlab -nodisplay -r "roi_tr_pattern_regression_SL_permPhase($ei)"
# matlab -nodisplay -r "roi_tr_temporal_regression_SLeach_betaClassification($ei)"

matlab -nodisplay -r "wholeBrain_tr_temporal_regression_SL($ei)"

# matlab -nodisplay -r "roi_tr_pattern_regression_LL_permPhase($perm)"


# matlab -nodisplay -r "roi_tr_bined_pattern_granger_SL_permSubj($s)"
# matlab -nodisplay -r "roi_tr_bined_pattern_regression_LL_permSubj($perm)"

# sqeueu -u huichuan
# scancel <jobid>









