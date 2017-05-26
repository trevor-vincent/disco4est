#!/bin/bash

function write_submit {
    cat <<EOF1 > submit.sh
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -o ${1}/disco4est.stdout
#PBS -e ${1}/disco4est.stderr
#PBS -d .
#PBS -S /bin/bash
#PBS -N $3

module purge;
module load gcc/4.8.1;
module load openmpi/gcc/1.8.3;
module load openblas/1.13-singlethreaded;
cd ${1}

time mpirun -np 8 ./${2}  2>&1 | tee disco4est.out
EOF1

}

function write_options {

    cat <<EOF > options.input
[p4est]
min_quadrants = -1
min_level = 0
fill_uniform = 1

[amr]
num_unifrefs = $1
num_randrefs = 0

[flux]
ip_flux_penalty = $2

[problem]
deg_R2 = $3
deg_quad_R2 = $4
deg_R1 = $5
deg_quad_R1 = $6
deg_R0 = $7
deg_quad_R0 = $7

R0 = 10
R1 = $8
R2 = $9

[solver]
ksp_type = gmres
ksp_atol = 1e-20
ksp_rtol = 1e-10
ksp_maxit = 10000
ksp_view = 1
ksp_monitor = 1
EOF
}

arr1=(0 1 2 ) #num_unifrefs
arr2=( 2 100 ) #penalty flux
arr3=( 8 15 ) #degR2
arr4=( 8 15 ) #degR2integ
arr5=( 4 8 ) #degR1
arr6=( 4 8 ) #degR1integ
arr7=( 2 5 ) #degR0
arr8=( 50 ) #R1
arr9=( 1000 10000000 10000000000 ) #R2

for a in "${arr1[@]}"
do
    for b in "${arr2[@]}"
    do
	for c in "${arr3[@]}"
	do
	    for d in "${arr4[@]}"
	    do
		for e in "${arr5[@]}"
		do
		    for f in "${arr6[@]}"
		    do
			for g in "${arr7[@]}"
			do
			    for h in "${arr8[@]}"
			    do
				for i in "${arr9[@]}"
				do
				
			    # for h in "${arr8[@]}"
			    # do
				#break for loop on Gauss offset if were not using Gauss
				# if [ "$g" -eq "0" ]; then
				# if [ "$h" -gt "0" ]; then
				#     break
				# fi
				# fi
 				NEWDIR="sphere_unif_${a}_pen_${b}_degR2_${c}_degR2int_${d}_degR1_${e}_degR1int_${f}__degR0_${g}_R1_${h}_R2_${i}_Solver_gmres"
				mkdir $NEWDIR
				cd $NEWDIR
				SHORTNAME="sph${a}${b}${c}${d}${e}${f}${g}${h}${i}"
				rundir=$PWD
				executable_path=$1
				executable=$2
				#	executable="constantdensitystar_driver"
				#	executable_path=""
				write_options $a $b $c $d $e $f $g $h $i
				write_submit $rundir $executable $SHORTNAME
				ln -s "${executable_path}/${executable}" "${PWD}/${executable}"
				#	qsub submit.sh
				cd ..
				done
			    done
			done
		    done  
		done
	    done
	done
    done  
done

