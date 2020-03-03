#!/bin/bash
#SBATCH -J 'roi2'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 3:00:00
#SBATCH -c 10
#SBATCH --array=1-61
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END
#SBATCH --mail-user=hcchang73@gmail.com


module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

perm=$SLURM_ARRAY_TASK_ID
declare -i perm


matlab -nodisplay -r "roi2rois_temporal_lagcorr_LL_leave1out_permPhase($perm)"










