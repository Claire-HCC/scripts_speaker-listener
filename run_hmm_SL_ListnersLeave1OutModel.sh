#!/bin/bash
#SBATCH -J 'hmm_ListnersLeave1OutModel'
#SBATCH -o slurm-%j.out
#SBATCH --time=1-0 # days-hours
#SBATCH -c 20 # 1 core per task
#SBATCH -n 1 # tasks
#SBATCH --nodes 2 # minimum of minnodes nodes
#SBATCH --mem-per-cpu=500M
# e.g. sbatch run_hmm_SL_ListnersLeave1OutModel.sh

module load pyger

#exps=("pieman" "bronx" "merlin" "sherlock")
expdir="/mnt/sink/scratch/claire/speaker-listener/"
exp="merlin"
fs=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat`
rname="aCUN"

python ./hmm_SL_ListnersLeave1OutModel.py $rname $exp

seff $SLURM_JOBID
# for f in $fs; do
    # rname=`echo ${i%.*} | cut -d'_' -f 3`
    # python ./hmm_SL_findOptinalEventN.py $rname $exp
# done

sacct -o reqmem,maxrss,averss,elapsed -j %j
# scontrol show config
# squeue -u huichuan





