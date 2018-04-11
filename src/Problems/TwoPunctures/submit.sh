#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -o /scratch/p/pfeiffer/tvincent/d4est/d4est_build/TwoPunctures/disco4est.stdout
#PBS -e /scratch/p/pfeiffer/tvincent/d4est/d4est_build/TwoPunctures/disco4est.stderr
#PBS -d .
#PBS -S /bin/bash
#PBS -N two_punctures

source /home/p/pfeiffer/tvincent/d4est.env
cd /scratch/p/pfeiffer/tvincent/d4est/d4est_build/TwoPunctures

time mpirun -np 8 ./two_punctures_driver  2>&1 | tee disco4est.out
