#!/usr/local/bin/bash -l

export LOCAL_RANK=$SLURM_LOCALID
export GLOBAL_RANK=$SLURM_PROCID
# The first 4 ranks get assigned to GPUs
# Subsequent ranks get no GPU
export GPUS=(3 2 1 0)
export NUMA=(0 1 2 3 0 1 2 3)
export NUMA_NODE=${NUMA[$LOCAL_RANK]}

export CUDA_VISIBLE_DEVICES=${GPUS[$SLURM_LOCALID%8]}


numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE bash -c "$@"

