#!/bin/bash
#SBATCH -J 'permPahse'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 5:00:00
#SBATCH -c 10
#SBATCH --array=1-1000
#SBATCH --mem-per-cpu=500M
#SBATCH --mail-type=END
#SBATCH --mail-user=hcchang73@gmail.com

module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

perm=$SLURM_ARRAY_TASK_ID
declare -i perm


matlab -nodisplay -r "roi_pattern_regression_SL_each($perm)"











