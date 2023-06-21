#! /bin/bash
set +eu

nvsmi_logger_PID=0
function kill_nvsmi()
{
    set +x
    if (( nvsmi_logger_PID != 0 ))
    then
        kill $nvsmi_logger_PID
    fi
}
trap kill_nvsmi ERR
trap kill_nvsmi EXIT

lrank=$OMPI_COMM_WORLD_RANK
mpi_total_procs=$OMPI_COMM_WORLD_SIZE
# need to check in run script that the variables make sense and are
# exported!
# (( compute_tasks = mpi_total_procs - io_tasks ))
compute_tasks=$mpi_total_procs

echo compute_tasks: $compute_tasks, lrank: $lrank

# Local Task 0 runs always the nvidia-smi logger 
if [[ "$lrank" == 0 ]] && [[ ${ENABLE_NVIDIA_SMI_LOGGER:-"no"} == "yes" ]]
then
    set +x
    basedir=${basedir:=.}
    # Start logger in background. It will be killed by the ERR trap or kill_nvsmi.
    loop_repetition_time=500000000 # in nano seconds
    while sleep 0.$(( ( 1999999999 - 1$(date +%N) ) % loop_repetition_time ))
    do
        LC_TIME=en_US date -Ins
        nvidia-smi --format=csv --query-gpu=index,power.draw,utilization.gpu,temperature.gpu,memory.used
    done > nvsmi.log.${lrank} &
    nvsmi_logger_PID=$!
    set -x
fi

# numanode=(3 1 7 5)
# nics=(mlx5_4 mlx5_5 mlx5_0 mlx5_1 mlx5_10 mlx5_11 mlx5_7 mlx5_8)
# reorder=(0 1 2 3 4 5 6 7)

gpus=(0 1 2 3 4 5 6 7)

# splitted the blocks for later use of gpnl nodes as I/O server
if (( lrank < compute_tasks ))
then
    echo Compute process $OMPI_COMM_WORLD_RANK on $(hostname)
    export CUDA_VISIBLE_DEVICES=${gpus[$((lrank % ${#gpus[@]} ))]}
else
    echo IO process $OMPI_COMM_WORLD_RANK on $(hostname)
fi

export KMP_AFFINITY=scatter

$@
return=$?

echo "nvsmi_logger_PID $nvsmi_logger_PID"
kill_nvsmi
exit $return

