#!/bin/bash
#SBATCH -J 'matlab_job'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 03:00:00
#SBATCH -c 5
#SBATCH --array=1-2
#SBATCH --mem-per-cpu=1G


module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

lagi=$SLURM_ARRAY_TASK_ID
declare -i lagi

matlab -nodisplay -r "roi_tr_bined_pattern_granger_SL($lagi)"

# matlab -nodisplay -r "roi_tr_bined_pattern_regression_LL"



# matlab ./roi_tr_bined_pattern_granger_SL.m $lagi






