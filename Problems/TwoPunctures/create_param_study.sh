#!/bin/bash

function write_submit {
    cat <<EOF1 > submit.sh
#PBS -l nodes=${6}:ppn=8
#PBS -l walltime=$5:00:00
#PBS -o ${1}/disco4est.stdout
#PBS -e ${1}/disco4est.stderr
#PBS -d .
#PBS -S /bin/bash
#PBS -N $3

module purge;
source /home/p/pfeiffer/tvincent/d4est.env
cd ${1}

time mpirun -np $4 ./${2}  2>&1 | tee disco4est.out
EOF1

}

function write_options {

    cat <<EOF > options.input
[initial_grid]
min_quadrants = -1
min_level = $1
fill_uniform = 1
deg = $2
deg_quad = $2

[amr]
scheme = smooth_pred
num_of_amr_steps = 15
max_degree = 7 
gamma_h = $3
gamma_p = 0.1
gamma_n = 1.
percentile = $4

[flux]
name = sipg
sipg_penalty_prefactor = 2.0
sipg_flux_h = $5
sipg_penalty_fcn = maxp_sqr_over_minh


[geometry]
name = cubed_sphere
R0 = 10
R1 = 20
R2 = $6
compactify_outer_shell = $7
compactify_inner_shell = 0
DX_compute_method = analytic
JAC_compute_method = numerical


[d4est_vtk_geometry]
name = cubed_sphere
R0 = 1
R1 = 2
R2 = 3
compactify_outer_shell = 0
compactify_inner_shell = 0
DX_compute_method = analytic
JAC_compute_method = numerical

[problem]
do_not_solve = 0
deg_vol_quad_inc_inner = $8
deg_vol_quad_inc_outer = $9

[d4est_solver_newton]
imin = 3
imax = 3
monitor = 1
rtol = 1e-5
atol = 1e-15

[krylov_petsc]
ksp_type = fcg
ksp_atol = 1e-15
ksp_rtol = 1e-5
ksp_max_it = 10000
ksp_view = 1
ksp_monitor = 1
ksp_converged_reason = 1
ksp_initial_guess_nonzero = 0
ksp_monitor_singular_value = 1

[quadrature]
name = legendre
EOF
}

arr1=( 3 ) #initial degree
arr2=( 1 2 ) #percentile
arr3=(.25 2.0 6.0) #gammah
arr4=( 5 ) #penalty
arr5=( "H_EQ_J_DIV_SJ_MIN_LOBATTO" "H_EQ_FACE_DIAM" "H_EQ_J_DIV_SJ_QUAD") #hrefine til inview
arr6=( 40 100 1000 ) #domain size
arr7=( 0 ) #Gauss offset
arr8=( 4 )
arr9=( 4 )

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

 				NEWDIR="2pun_${a}_${b}_${c}_${d}_${e}_${f}_${g}_${h}_${i}"
				mkdir $NEWDIR
				cd $NEWDIR
				SHORTNAME="2pun${a}${b}${c}${d}${e}${f}${g}${h}${i}"
				rundir=$PWD
				executable_path=$1
				executable=$2
				cores=$3
				hours=$4
				nodes=$((${cores} / 8))
				write_options $a $b $c $d $e $f $g $h $i
				write_submit $rundir $executable $SHORTNAME $cores $hours $nodes
				ln -s "${executable_path}/${executable}" "${PWD}/${executable}"
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
