#!/bin/bash
#SBATCH -J 'matlab_roi'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 1:00:00
#SBATCH -c 10
#SBATCH --array=1-61
#SBATCH --mem-per-cpu=500M


module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

perm=$SLURM_ARRAY_TASK_ID
declare -i perm


matlab -nodisplay -r "roi_tr_temporal_lagcorr_SLg_permPhase($perm)"










