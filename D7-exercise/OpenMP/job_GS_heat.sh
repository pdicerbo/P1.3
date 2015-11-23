#!/bin/bash
#PBS -N MHPC_OPENMP_HEAT_GS
#PBS -l walltime=02:00:00 
#PBS -l nodes=1:ppn=20:gpu  
#PBS -T flush_cache 
#PBS -q reserved4

module purge
module load intel

cd $PBS_O_WORKDIR

N_STEPS=100
DT=0.0001

echo -e "\n\n\n\t==== Start Serial Execution ===="
./GS_heat_openmp_serial.x ${N_STEPS} ${DT}
echo -e "\n\t==== End Serial Execution ===="

echo -e "\n\n\n\t==== Start Parallel Executions ===="
for THREADS in 2 4 8 16;
do
 export OMP_NUM_THREADS=$THREADS
 ./GS_heat_openmp.x ${N_STEPS} ${DT}
done
echo -e "\n\t==== End Parallel Execution ====\n\n"

