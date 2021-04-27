#!/bin/bash
#SBATCH -J 'matlab_job'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 5:00:00
#SBATCH -c 10
#SBATCH --array=1-61
#SBATCH --mem-per-cpu=500M
#SBATCH --mail-type=END
#SBATCH --mail-user=hcchang73@gmail.com

module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

perm=$SLURM_ARRAY_TASK_ID
declare -i perm


matlab -nodisplay -r "roi2rois_temporal_circularlagcorr_LL_selfself($perm)"












