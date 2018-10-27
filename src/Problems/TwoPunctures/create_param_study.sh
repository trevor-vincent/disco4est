#!/bin/bash

function write_submit_scinet {
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

function write_submit_graham {
    cat <<EOF1 > submit.sh
#!/bin/bash
#SBATCH --ntasks=${4}               # number of MPI processes
#SBATCH --mem-per-cpu=1024M      # memory; default unit is megabytes
#SBATCH --time=$5-00:00           # time (DD-HH:MM)

source /home/tvincent/d4est.env
cd ${1}
time mpirun -np $4 ./${2}  2>&1 | tee disco4est.out

EOF1
}


function write_submit_minerva {
    cat <<EOF1 > submit.sh

#!/bin/bash -
#SBATCH -J testd4est           # Job Name
#SBATCH -o testd4est.stdout    # Output file name
#SBATCH -e testd4est.stderr    # Error file name
#SBATCH -n 1                    # Number of cores
#SBATCH --ntasks-per-node 16    # DO NOT TOUCH: number of MPI ranks per node
#SBATCH -t 24:0:00               # Max run time
#SBATCH --no-requeue

source path/to/Support/Scripts/submit.env
mpirun -n 4 path/to/problem_executable

EOF1
}

function write_submit_niagara {
    cat <<EOF1 > submit.sh
#!/bin/bash
#SBATCH --nodes=${4}
#SBATCH --time=$5-00:00           # time (DD-HH:MM)

source /home/p/pfeiffer/tvincent/d4est.env
cd ${1}
time mpirun -np $4 --map-by node ./${2}  2>&1 | tee disco4est.out

EOF1
}


function write_options {

    cat <<EOF > options.input
[initial_mesh]
min_quadrants = -1
min_level = 2
fill_uniform = 1
region0_deg = $1
region0_deg_quad_inc = $2
region1_deg = $1
region1_deg_quad_inc = $2
region2_deg = $1
region2_deg_quad_inc = $2
load_from_checkpoint = 0
load_from_newton_checkpoint = 1
newton_checkpoint_prefix = checkpoint_newton_2
checkpoint_type = D4EST_CHKPT_HISTORY_H5
checkpoint_prefix = checkpoint
checkpoint_number = 6
initial_checkpoint_number = 1

[mesh_parameters]
face_h_type = FACE_H_EQ_J_DIV_SJ_QUAD
volume_h_type = VOL_H_EQ_DIAM
max_degree = 20

[problem]
do_not_solve = 0
use_compactified_size_params = 0
use_error_l2_as_estimator = 0
use_dirichlet = 1
use_new_multigrid_level_guesser = 0

[amr]
scheme = smooth_pred
num_of_amr_steps = 100
max_degree = 10
gamma_h = $3
gamma_p = 0.1
gamma_n = 1.
inflation_size = 128
percentile = $4
initial_pred = 10000

[flux]
name = sipg
sipg_penalty_prefactor = $5
sipg_penalty_fcn = maxp_sqr_over_minh

[geometry]
name = cubed_sphere
R0 = 10
R1 = 20
R2 = $6
compactify_outer_shell = 1
compactify_inner_shell = 0
DX_compute_method = analytic
JAC_compute_method = numerical

[compactified_geometry]
name = cubed_sphere
R0 = 1
R1 = 2
R2 = 3
compactify_outer_shell = 1
compactify_inner_shell = 0
DX_compute_method = analytic
JAC_compute_method = numerical

[d4est_solver_newton]
atol = 1e-15
rtol = 1e-15
imax = 1000000
imin = 1
monitor = 1

[d4est_solver_fcg]
atol = 1e-15
rtol = 1e-15
iter = 1000000
monitor = 1
vi = 5
precond_flag = 0

[d4est_solver_cg]
atol = 1e-15
rtol = 1e-15
iter = 1000000000
monitor = 1
[d4est_vtk]
filename = two_punctures
geometry_section = geometry
output_type = ascii
grid_type = dg
write_tree = 1
write_level = 1
write_rank = 1
wrap_rank = 0
write_deg = 1

[d4est_vtk_compactified]
filename = two_punctures_compactified
geometry_section = compactified_geometry
output_type = ascii
grid_type = dg
write_tree = 1
write_level = 1
write_rank = 1
wrap_rank = 0
write_deg = 1

[d4est_vtk_compactified_corner]
filename = two_punctures_compactified_corner
geometry_section = compactified_geometry
output_type = ascii
grid_type = corner
write_tree = 1
write_level = 1
write_rank = 1
wrap_rank = 0
write_deg = 1

[d4est_vtk_corner]
filename = two_punctures_corner
geometry_section = geometry
output_type = ascii
grid_type = corner
write_tree = 1
write_level = 1
write_rank = 1
wrap_rank = 0
write_deg = 1

[quadrature]
name = legendre

[d4est_solver_newton_petsc]
snes_atol = 1e-15
snes_rtol = 1e-50
snes_stol = 1e-15
snes_max_funcs = 1000000000
snes_type = newtonls
snes_max_it = 5
snes_monitor = 1
snes_linesearch_order = 3
snes_linesearch_monitor = 1
snes_converged_reason = 1
snes_view = 1
checkpoint_every_n_newton_its = 1

[d4est_solver_krylov_petsc]
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
vcycle_imax = $7;
vcycle_rtol = 1e-9;
vcycle_atol = 0.;
smoother_name = mg_smoother_cheby
bottom_solver_name = $8
use_analyzer = 0

[mg_bottom_solver_cg]
bottom_iter = 100;
bottom_rtol = 1e-10;
bottom_atol = 0.;
bottom_print_residual_norm = 0;

[mg_smoother_cheby]
cheby_imax = $9;
cheby_eigs_cg_imax = 30;
cheby_eigs_lmax_lmin_ratio = 30.;
cheby_eigs_max_multiplier = 1.1;
cheby_eigs_reuse_fromdownvcycle = 0;
cheby_eigs_reuse_fromlastvcycle = 0;
cheby_print_residual_norm = 0;
cheby_print_eigs = 0;
cheby_use_new_cg_eigs = 1;

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

arr1=( 1 2 ) #min_level
arr2=( 0 3 ) #deg
arr3=( .25 .75 1 2 5) #deg_quad_inc
arr4=( 10 25 ) #hrefine til inview
arr5=( 2 100 ) #penalty
arr6=( 10000000000 100000000000 ) #domain size
arr7=( 1 5 ) #Gauss offset
arr8=( "mg_smoother_cheby" "mg_bottom_solver_cg" )
arr9=( 15 )

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
				write_submit_niagara $rundir $executable $SHORTNAME $cores $hours $nodes
				cp "${executable_path}/${executable}" "${PWD}/${executable}"
				cp "${executable_path}/*.conf" "${PWD}/"
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

