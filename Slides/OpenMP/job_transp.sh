#!/bin/bash
#PBS -N MHPC_OPENMP_TRASP
#PBS -l walltime=02:00:00 
#PBS -l nodes=1:ppn=20:gpu  
#PBS -T flush_cache 
#PBS -q reserved4

module purge
module load intel

cd $PBS_O_WORKDIR

M_SIZE=16384

echo -e "\n\n\n\t==== Start Serial Execution ===="
./fast_transpose_openmp_serial.x ${M_SIZE}
echo -e "\n\t==== End Serial Execution ===="

echo -e "\n\n\n\t==== Start Parallel Executions ===="
for THREADS in 2 4 8 16;
do
 export OMP_NUM_THREADS=$THREADS
 ./fast_transpose_openmp.x ${M_SIZE}
done
echo -e "\n\t==== End Parallel Execution ====\n\n"

