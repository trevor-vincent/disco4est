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

function write_submit_cedar {
    cat <<EOF1 > submit.sh
#!/bin/bash
#SBATCH --ntasks=${4}               # number of MPI processes
#SBATCH --time=$5-00:00           # time (DD-HH:MM)

source /home/tvincent/d4est.env
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
#SBATCH -J test_d4est           # Job Name
#SBATCH -o test_d4est.stdout    # Output file name
#SBATCH -e test_d4est.stderr    # Error file name
#SBATCH -n 1                    # Number of cores
#SBATCH --ntasks-per-node 16    # DO NOT TOUCH: number of MPI ranks per node
#SBATCH -t 24:0:00               # Max run time
#SBATCH --no-requeue

source path/to/Support/Scripts/submit.env
mpirun -n 4 path/to/problem_executable

EOF1
}


function write_options {

    cat <<EOF > options.input
[initial_mesh]
min_quadrants = -1
min_level = $1
fill_uniform = 1
region0_deg = 1
region0_deg_quad_inc = $2
region1_deg = 1
region1_deg_quad_inc = 0

[problem]
use_dirichlet = $3

[flux]
name = sipg
sipg_penalty_prefactor = $4
sipg_flux_h = $5
;H_EQ_VOLUME_DIV_AREA
sipg_penalty_fcn = $6
;maxp_sqr_over_minh

[amr]
scheme = uniform_p
num_of_amr_steps = 200
max_degree = 7


[geometry]
name = cubed_sphere_7tree
R0 = 10
R1 = $7
compactify_outer_shell = 0
compactify_inner_shell = $8
DX_compute_method = analytic
JAC_compute_method = numerical

[d4est_vtk]
filename = poisson_exp1ox_uniform
geometry_section = geometry
output_type = ascii
grid_type = dg
write_tree = 1
write_level = 1
write_rank = 1
wrap_rank = 0
write_deg = 1

[d4est_solver_cg]
iter = 1000000
atol = 1e-15
rtol = 1e-15
monitor = 1

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

