#!/bin/bash
#SBATCH -J 'Exp'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 2:00:00
#SBATCH -c 10
#SBATCH --array=1-4
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END
#SBATCH --mail-user=hcchang73@gmail.com

module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

ei=$SLURM_ARRAY_TASK_ID
declare -i ei


matlab -nodisplay -r "roi_pattern_lagcorr_LL_leave1out($ei)"


# sqeueu -u huichuan
# scancel <jobid>









