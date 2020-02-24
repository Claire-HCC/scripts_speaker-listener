#!/bin/bash
#SBATCH -J 'matlab_job'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 3:00:00
#SBATCH -c 10
#SBATCH --array=11-18
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=NONE
#SBATCH --mail-user=hcchang73@gmail.com

module load fsl

expdir="/mnt/sink/scratch/claire/speaker-listener/"

# fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

perm=$SLURM_ARRAY_TASK_ID
declare -i perm

flirt -in /scratch/claire/speaker-listener/merlin/fmri/timeseries/tr/wholeBrain/3x3x3mm/listener$(printf "%02d" "$perm").nii -ref /scratch/claire/speaker-listener/roi_mask/MNI152NLin2009cAsym_3x3x4mm_brain.nii -out /scratch/claire/speaker-listener/merlin/fmri/timeseries/tr/wholeBrain/listener$(printf "%02d" "$perm").nii  -applyxfm -usesqform
fslchfiletype NIFTI /scratch/claire/speaker-listener/merlin/fmri/timeseries/tr/wholeBrain/listener$(printf "%02d" "$perm").nii.gz /scratch/claire/speaker-listener/merlin/fmri/timeseries/tr/wholeBrain/listener$(printf "%02d" "$perm").nii













