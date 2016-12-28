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
min_level = 1
fill_uniform = 1

[amr]
degmax = $9
amr_levels = 30
amr_inflation_size = $8
initial_degree = $1
percentile = $2
gamma_h = $3
gamma_p = 0.1

[flux]
ip_flux_penalty = $4

[problem]
hrefine_til_inview = $5
domain_size = $6 
use_Gauss_integ = 1
deg_offset_for_Gauss_integ = $7

[solver]
krylov_type = cg
EOF
}

arr1=( 1 2 ) #initial degree
arr2=( 5 10 ) #percentile
arr3=(.25 2.0 6.0) #gammah
arr4=( 2 ) #penalty
arr5=( 1 ) #hrefine til inview
arr6=( 20000 1000 10 ) #domain size
arr7=( 1 2 ) #Gauss offset
arr8=( 64 128 256 )
arr9=( 7 5 )

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
 				NEWDIR="2pun_deg_${a}_perc_${b}_gamh_${c}_pen_${d}_hreftilinview_${e}_domsize_${f}_useGauss_1_Gaussoffset_${g}_AmrInfSize_${h}_Solver_cg_degmax_${i}"
				mkdir $NEWDIR
				cd $NEWDIR
				SHORTNAME="2pun${a}${b}${c}${d}${e}${f}${g}${h}${i}"
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

