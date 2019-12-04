#!/bin/bash
#SBATCH -J 'hmm_SL_findOptinalEventN'
#SBATCH -o slurm-%j.out
#SBATCH --time=1-0 # days-hours
#SBATCH -c 2 # 1 core per task
#SBATCH -n 2 # tasks
#SBATCH --nodes 2 # minimum of minnodes nodes
#SBATCH --mem-per-cpu=4G
# this takes 2,5 hr/roi

module load pyger

#exps=("pieman" "bronx" "merlin" "sherlock")
expdir="/mnt/sink/scratch/claire/speaker-listener/"
exp="bronx"
fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`

for f in $fs; do
    rname=`echo ${i%.*} | cut -d'_' -f 3`
    python ./hmm_SL_findOptinalEventN.py $rname $exp
done
# python ./hmm_SL_findOptinalEventN.py dPCC
# sacct -o reqmem,maxrss,averss,elapsed -j 13825317
# scontrol show config
# jupyter nbconvert --to python hmm_SL_findOptinalEventN.ipynb
# squeue -u huichuan





