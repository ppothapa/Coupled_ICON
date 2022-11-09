#!/usr/local/bin/bash -l

export LOCAL_RANK=$SLURM_LOCALID
export GLOBAL_RANK=$SLURM_PROCID
export NUMA_NODE=$(($LOCAL_RANK/16))

numactl --physcpubind=$LOCAL_RANK --membind=$NUMA_NODE bash -c "$@"