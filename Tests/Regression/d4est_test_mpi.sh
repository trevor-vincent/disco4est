echo "Starting MPI Regression test"
echo ${PWD}
sed 's|num_of_amr_steps = .*|num_of_amr_steps = 0|g' ../ConstantDensityStar/options.input > mpi_options.input
mpirun -np 4 ../ConstantDensityStar/constant_density_star_driver mpi_options.input > disco4est.out
RES0=$(cat disco4est.out | grep -c "0.0000096078")
RES1=$(cat disco4est.out | grep -c "Completed problem")
RES=$(echo $RES0 + $RES1 | bc) 
rm mpi_options.input
if [ $RES -ne "2" ]; then
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
    rm disco4est.out
fi



