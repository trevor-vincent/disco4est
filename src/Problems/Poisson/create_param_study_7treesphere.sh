#!/bin/bash

function write_submit_graham {
    cat <<EOF1 > submit.sh
#!/bin/bash
#SBATCH --ntasks=${4}               # number of MPI processes
#SBATCH --mem-per-cpu=1024M      # memory; default unit is megabytes
#SBATCH --time=$5-00:00           # time (DD-HH:MM)

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK 
source /home/tvincent/d4est.env
cd ${1}
time mpirun -np $4 ./${2}  2>&1 | tee disco4est.out

EOF1
}

function write_options {

    cat <<EOF > options.input
[initial_mesh]
min_quadrants = -1
min_level = 0
fill_uniform = 1
region0_deg = $1
region0_deg_quad_inc = $2
region1_deg = $1
region1_deg_quad_inc = $2

[problem]
use_dirichlet = 1

[flux]
name = sipg
sipg_penalty_prefactor = $3
sipg_flux_h = $4
sipg_penalty_fcn = maxp_sqr_over_minh

[amr]
scheme = uniform_h
num_of_amr_steps = 100
max_degree = 7

[geometry]
name = cubed_sphere_7tree
R0 = $5
R1 = $6
compactify_outer_shell = 0
compactify_inner_shell = 0
DX_compute_method = analytic
JAC_compute_method = numerical

[d4est_vtk]
filename = poisson_sinx_uniform
geometry_section = geometry
output_type = ascii
grid_type = dg
write_tree = 1
write_level = 1
write_rank = 1
wrap_rank = 0
write_deg = 1

[quadrature]
name = legendre

[d4est_solver_krylov_petsc]
ksp_type = fcg
ksp_atol = 5e-15
ksp_rtol = 1e-20
ksp_max_it = 10000
ksp_view = 0
ksp_monitor = 1
ksp_converged_reason = 1
ksp_initial_guess_nonzero = 0
ksp_monitor_singular_value = 0

[multigrid]
vcycle_imax = $7;
vcycle_rtol = $8;
vcycle_atol = 0.;
smoother_name = mg_smoother_cheby
bottom_solver_name = mg_bottom_solver_cg

[mg_bottom_solver_cg]
bottom_iter = $9;
bottom_rtol = 1e-10;
bottom_atol = 0.;
bottom_print_residual_norm = 0;

[mg_smoother_cheby]
cheby_imax = 8;
cheby_eigs_cg_imax = 30;
cheby_eigs_lmax_lmin_ratio = 30.;
cheby_eigs_max_multiplier = 1.;
cheby_eigs_reuse_fromdownvcycle = 0;
cheby_eigs_reuse_fromlastvcycle = 0;
cheby_print_residual_norm = 0;
cheby_print_eigs = 0;

[mg_bottom_solver_cheby]
cheby_imax = 15;
cheby_eigs_cg_imax = 30;
cheby_eigs_lmax_lmin_ratio = 30.;
cheby_eigs_max_multiplier = 1.;
cheby_eigs_reuse_fromdownvcycle = 0;
cheby_eigs_reuse_fromlastvcycle = 0;
cheby_print_residual_norm = 0;
cheby_print_eig = 0;


EOF
}

arr1=( 1 2 3 4 5 ) 
arr2=( 0 1 ) #percentile
arr3=( 2. 10. ) #gammah
arr4=( "H_EQ_VOLUME_DIV_AREA" "H_EQ_J_DIV_SJ_MIN_LOBATTO" ) #penalty
arr5=( .5 ) #hrefine til inview
arr6=( 1 ) #domain size
arr7=( 1 3 ) #Gauss offset
arr8=( 1e-15 )
arr9=( 100 )

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

