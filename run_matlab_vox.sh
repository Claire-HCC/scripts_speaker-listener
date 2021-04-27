#!/bin/bash
#SBATCH -J 'matlab_job'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 4:00:00
#SBATCH -c 10
#SBATCH --array=2-17
#SBATCH --mem-per-cpu=500M
#SBATCH --mail-type=END
#SBATCH --mail-user=hcchang73@gmail.com

module load matlab

expdir="/mnt/sink/scratch/claire/speaker-listener/"

perm=$SLURM_ARRAY_TASK_ID
declare -i perm


matlab -nodisplay -r "vox2voxs_temporal_circularlagcorr_LL_gg($perm)"












