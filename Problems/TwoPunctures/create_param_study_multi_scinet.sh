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
[initial_mesh]
min_quadrants = -1
min_level = $1
fill_uniform = 1
region0_deg = $2
region0_deg_quad_inc = $3
region1_deg = $2
region1_deg_quad_inc = 0

[amr]
scheme = smooth_pred
num_of_amr_steps = 15
max_degree = 7 
gamma_h = .25
gamma_p = 0.1
gamma_n = 1.
inflation_size = 128
percentile = 5
amr_level_for_uniform_p = $4
use_puncture_finder = $5

[problem]
do_not_solve = 0
number_of_punctures = 3
puncture0_M = 0.269115643269342
puncture0_X = 0.015239564243350
puncture0_Y = -0.693307138836620
puncture0_Z = 0.000000000000000
puncture0_PX = 0.058489108578530
puncture0_PY = 0.008229897640753
puncture0_PZ = 0.000000000000000
puncture0_SX = 0.000000000000000
puncture0_SY = 0.000000000000000
puncture0_SZ = -0.013354678737630
puncture1_M = 0.406384234358315
puncture1_X = -2.316041324900483
puncture1_Y = 1.827477408027033
puncture1_Z = 0.000000000000000
puncture1_PX = -0.028489062296454
puncture1_PY = -0.149654960981410
puncture1_PZ = 0.000000000000000
puncture1_SX = 0.000000000000000
puncture1_SY = 0.000000000000000
puncture1_SZ = -0.133191216519657
puncture2_M = 0.324500122372343
puncture2_X = -1.027934858588472
puncture2_Y = -2.271160404789802
puncture2_Z = 0.000000000000000
puncture2_PX = 0.164072286507148
puncture2_PY = 0.051549165440513
puncture2_PZ = 0.000000000000000
puncture2_SX = 0.000000000000000
puncture2_SY = 0.000000000000000
puncture2_SZ = -0.170845803046062

[flux]
name = sipg
sipg_penalty_prefactor = $6
sipg_flux_h = H_EQ_J_DIV_SJ_MIN_LOBATTO
sipg_penalty_fcn = maxp_sqr_over_minh

[geometry]
name = cubed_sphere_7tree
R0 = $7
R1 = $8
compactify_outer_shell = 0
compactify_inner_shell = $9
DX_compute_method = analytic
JAC_compute_method = numerical


[d4est_vtk_geometry]
name = cubed_sphere_7tree
R0 = .5
R1 = 1
compactify_outer_shell = 0
compactify_inner_shell = $9
DX_compute_method = analytic
JAC_compute_method = numerical

[quadrature]
name = legendre

[newton_petsc]
snes_atol = 1e-15
snes_rtol = 1e-50
snes_stol = 1e-1
snes_max_funcs = 1000000000
snes_type = newtonls
snes_max_it = 5
snes_monitor = 1
snes_linesearch_order = 3
snes_linesearch_monitor = 1
snes_converged_reason = 1
snes_view = 1
snes_ksp_ew = 1

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

[multigrid]
vcycle_imax = 1;
vcycle_rtol = 1e-9;
vcycle_atol = 0.;
smoother_name = mg_smoother_cheby
bottom_solver_name = mg_bottom_solver_cheby

[mg_bottom_solver_cg_d4est]
bottom_iter = 100;
bottom_rtol = 1e-10;
bottom_atol = 0.;
bottom_print_residual_norm = 0;

[mg_smoother_cheby]
cheby_imax = 15;
cheby_eigs_cg_imax = 15;
cheby_eigs_lmax_lmin_ratio = 30.;
cheby_eigs_max_multiplier = 1.;
cheby_eigs_reuse_fromdownvcycle = 0;
cheby_eigs_reuse_fromlastvcycle = 0;
cheby_print_residual_norm = 0;
cheby_print_eigs = 0;

[mg_bottom_solver_cheby]
cheby_imax = 15;
cheby_eigs_cg_imax = 15;
cheby_eigs_lmax_lmin_ratio = 30.;
cheby_eigs_max_multiplier = 1.;
cheby_eigs_reuse_fromdownvcycle = 0;
cheby_eigs_reuse_fromlastvcycle = 0;
cheby_print_residual_norm = 0;
cheby_print_eig = 0;



EOF
}

arr1=( 2 3 ) #min_level
arr2=( 1 2 ) #deg
arr3=( 0 2) #deg_quad_inc
arr4=( 4 5 6 ) #hrefine til inview
arr5=( 0 1 ) #penalty
arr6=( 2.0 20.0 100.0 ) #domain size
arr7=( 10 12 ) #Gauss offset
arr8=( 1000 )
arr9=( 0 1 )

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
				cp "${executable_path}/${executable}" "${PWD}/${executable}"
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

