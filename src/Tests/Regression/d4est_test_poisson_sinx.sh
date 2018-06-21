echo "Starting PoissonSinx Regression test"
echo ${PWD}
cp ../Poisson/options.input poisson_sinx_options.input
../Poisson/poisson_sinx_uniform_driver poisson_sinx_options.input > disco4est.out
RES=$(cat disco4est.out | grep -c "64 512 512 0.0426895934625")
rm poisson_sinx_options.input
if [ $RES -ne "1" ]; then
    echo "Test Fail"
    echo "****************** POISSON SINX REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** POISSON SINX REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** POISSON SINX REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** POISSON SINX REGRESSION TEST OUTPUT BEGIN *****************"
    cat disco4est.out
    echo "****************** POISSON SINX REGRESSION TEST OUTPUT END *****************"
    echo "****************** POISSON SINX REGRESSION TEST OUTPUT END *****************"
    echo "****************** POISSON SINX REGRESSION TEST OUTPUT END *****************"
    echo "****************** POISSON SINX REGRESSION TEST OUTPUT END *****************"
    rm disco4est.out
    exit 1
else
    echo "Poisson Sinx Regression Test Success"
    rm disco4est.out
fi



