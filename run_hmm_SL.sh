#!/bin/bash
#SBATCH -J 'hmm_ListnersLeave1OutModel'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 6:00:00
#SBATCH -c 10
#SBATCH --array=1-61
# e.g. sbatch run_hmm_SL_ListenersLeave1Out_withinAcross.sh

module load pyger

# exps=("pieman" "bronx" "merlin" "sherlock")
expdir="/mnt/sink/scratch/claire/speaker-listener/"
exp="bronx"
rnames=`ls $expdir/$exp/fmri/timeseries/tr/roi/mor/zscore_listenerAll_*.mat | cut -d"/" -f14 | grep -o '^[^\.]*' | sed 's/zscore_listenerAll_//'`
# rname="aANG_L"

ri=$SLURM_ARRAY_TASK_ID

echo "Array Allocation Number: $SLURM_ARRAY_JOB_ID"

python ./hmm_SL_ListenersLeave1Out_withinAcross.py "$exp" "$ri"
# python ./hmm_SL_findListenersEventInSpeaker.py "$exp" "$ri"

echo "done!"
# sacct -o reqmem,maxrss,averss,elapsed -j %j
# scontrol show config
# squeue -u huichuan





