#PBS -l nodes=8:ppn=8
#PBS -l walltime=4:00:00
#PBS -o ${1}/disco4est.stdout
#PBS -e ${1}/disco4est.stderr
#PBS -d .
#PBS -S /bin/bash
#PBS -N 2punsphere

module purge;
module load gcc/4.8.1;
module load openmpi/gcc/1.8.3;
module load openblas/1.13-singlethreaded;
cd ${1}

time mpirun -np 64 ./${2}  2>&1 | tee disco4est.out
