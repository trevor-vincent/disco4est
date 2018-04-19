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

[problem]
deg_quad_inc = 0
eval_method = EVAL_BNDRY_FCN_ON_LOBATTO

[amr]
scheme = uniform_p
num_of_amr_steps = 5
max_degree = 7

[flux]
name = sipg
sipg_penalty_prefactor = $3
sipg_flux_h = $4
sipg_penalty_fcn = maxp_sqr_over_minh

[geometry]
name = cubed_sphere
R0 = .33333333333333333333333
R1 = .66666666666666666666666
R2 = $5
compactify_outer_shell = 0
compactify_inner_shell = 0
DX_compute_method = analytic
JAC_compute_method = numerical

[d4est_vtk_geometry]
name = cubed_sphere
R0 = .33333333333333333333333
R1 = .66666666666666666666666
R2 = $5
compactify_outer_shell = 0
compactify_inner_shell = 0
DX_compute_method = analytic
JAC_compute_method = numerical

[krylov_petsc]
ksp_type = cg
ksp_atol = 1e-15
ksp_rtol = 1e-20
ksp_max_it = 1000000
ksp_view = 0
ksp_monitor = 1
ksp_converged_reason = 1
ksp_initial_guess_nonzero = 0
ksp_monitor_singular_value = 1

[quadrature]
name = legendre

EOF
}

arr1=( 0 1 2 3 4 ) 
arr2=( 1 ) #percentile
arr3=( 20.0) #gammah
arr4=( "H_EQ_VOLUME_DIV_AREA" "H_EQ_J_DIV_SJ_MIN_LOBATTO" ) #penalty
arr5=( 1 3) #hrefine til inview
arr6=( 0 ) #domain size
arr7=( 0 ) #Gauss offset
arr8=( 0 )
arr9=( 0 )

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
				SHORTNAME="test${a}${b}${c}${d}${e}${f}${g}${h}${i}"
				rundir=$PWD
				executable_path=$1
				executable=$2
				hours=$4
				nodes_temp=$(($a * 13))
				nodes=$(( $nodes_temp > 1 ? $nodes_temp : 1 ))
				cores=$(($nodes * 8))
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

