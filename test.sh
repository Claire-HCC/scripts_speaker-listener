#!/bin/bash
#Author: Anne
#Purpose: runpcit first on simulated data

#SBATCH -J 'pcit'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 7:00:00
#SBATCH -c 4
#SBATCH --array=41-42

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

echo "Array Allocation Number: $SLURM_ARRAY_JOB_ID"
echo "Array Index: $SLURM_ARRAY_TASK_ID"

module load matlab
boot_num=$SLURM_ARRAY_TASK_ID