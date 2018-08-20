echo "Starting MPI Regression test"
echo ${PWD}
sed 's|num_of_amr_steps = .*|num_of_amr_steps = 0|g' ../ConstantDensityStar/options.input > mpi_options.input
mpirun -np 8 ../ConstantDensityStar/constant_density_star_driver mpi_options.input > disco4est.out
RES0=$(cat disco4est.out | grep -c "0.0000096078")
RES1=$(cat disco4est.out | grep -c "Completed problem")
# mpirun -np 2 ../ConstantDensityStar/constant_density_star_driver mpi_options.input > disco4est.out
# RES2=$(cat disco4est.out | grep -c "0.0000096078")
mpirun -np 1 ../ConstantDensityStar/constant_density_star_driver mpi_options.input > disco4est.out
RES3=$(cat disco4est.out | grep -c "0.0000096078")

RES=$(echo $RES0 + $RES1 + $RES3 | bc) 
rm mpi_options.input
if [ $RES -ne "3" ]; then
    echo "Test Fail"
    echo "****************** MPI REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** MPI REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** MPI REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** MPI REGRESSION TEST OUTPUT BEGIN *****************"
    cat disco4est.out
    echo "****************** MPI REGRESSION TEST OUTPUT END *****************"
    echo "****************** MPI REGRESSION TEST OUTPUT END *****************"
    echo "****************** MPI REGRESSION TEST OUTPUT END *****************"
    echo "****************** MPI REGRESSION TEST OUTPUT END *****************"
    rm disco4est.out
    exit 1
else
    echo "MPI Regression Test Success"
    # echo "****************** MPI REGRESSION TEST OUTPUT BEGIN *****************"
    # echo "****************** MPI REGRESSION TEST OUTPUT BEGIN *****************"
    # echo "****************** MPI REGRESSION TEST OUTPUT BEGIN *****************"
    # echo "****************** MPI REGRESSION TEST OUTPUT BEGIN *****************"
    # cat disco4est.out
    # echo "****************** MPI REGRESSION TEST OUTPUT END *****************"
    # echo "****************** MPI REGRESSION TEST OUTPUT END *****************"
    # echo "****************** MPI REGRESSION TEST OUTPUT END *****************"
    # echo "****************** MPI REGRESSION TEST OUTPUT END *****************"
    rm disco4est.out
fi



