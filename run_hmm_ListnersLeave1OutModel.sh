#!/bin/bash
#SBATCH -J 'hmm_ListnersLeave1OutModel'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 2:00:00
#SBATCH -c 5
#SBATCH --array=10-120
# This takes around 12hr for 61 Mor ROIs, 111 Ks, for Merlin
# e.g. sbatch run_hmm_ListnersLeave1OutModel.sh

module load pyger

# exps=("pieman" "bronx" "merlin" "sherlock")
expdir="/scratch/claire/speaker-listener/"
exp="merlin"
rnames=`ls $expdir/$exp/fmri/timeseries/tr/network/mor/zscore_listenerAll_iscmasked_*.mat | cut -d"/" -f12 | grep -o '^[^\.]*' | sed 's/zscore_listenerAll_iscmasked_//'`
# rname="aANG_L"

K=$SLURM_ARRAY_TASK_ID

echo "Array Allocation Number: $SLURM_ARRAY_JOB_ID"

for rname in $rnames; do
    python ./hmm_ListnersLeave1OutModel.py "$exp" "$rname" $K
done

echo "done!"
# sacct -o reqmem,maxrss,averss,elapsed -j %j
# scontrol show config
# squeue -u huichuan





