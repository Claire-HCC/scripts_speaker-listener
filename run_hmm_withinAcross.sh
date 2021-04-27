#!/bin/bash
#SBATCH -J 'hmm_ListnersLeave1OutModel'
#SBATCH -o slurm-%j.out
#SBATCH -p all
#SBATCH -t 2:00:00
#SBATCH -c 5
#SBATCH --array=1-6
# e.g. sbatch run_hmm_SL_ListenersLeave1Out_withinAcross.sh

module load pyger

# exps=("pieman" "bronx" "merlin" "sherlock")
expdir="/scratch/claire/speaker-listener/"
exp="merlin"
rnames=`ls $expdir/$exp/fmri/timeseries/tr/network/mor/zscore_listenerAll_iscmasked_*.mat | cut -d"/" -f12 | grep -o '^[^\.]*' | sed 's/zscore_listenerAll_iscmasked_//'`
# rname="DMN2"

ri=$SLURM_ARRAY_TASK_ID

echo "Array Allocation Number: $SLURM_ARRAY_JOB_ID"

# python ./hmm_ListenersLeave1Out_withinAcross.py "$exp" "$ri"
# python ./hmm_ListenersLeave1Out_withinAcross_perm.py "$exp" "$ri"
# python ./hmm_findOptinalEventN.py 
python ./hmm_List.py "$exp" "$ri"

echo "done!"
# sacct -o reqmem,maxrss,averss,elapsed -j %j
# scontrol show config
# squeue -u huichuan





