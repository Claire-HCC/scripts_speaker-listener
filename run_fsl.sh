#!/bin/bash
#SBATCH -J 'matlab_job'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 5:00:00
#SBATCH -c 10
#SBATCH --array=1-48
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END
#SBATCH --mail-user=hcchang73@gmail.com

module load fsl

expdir="/mnt/sink/scratch/claire/speaker-listener/"


perm=$SLURM_ARRAY_TASK_ID
declare -i perm

dir_current=$(printf "/scratch/claire/speaker-listener/fmri_preprocesing/subjects/sub-%02d" "$perm")
echo "${dir_current}"
cd "${dir_current}"
# srun analyze.sh

srun ./scripts/apply-transform.sh









