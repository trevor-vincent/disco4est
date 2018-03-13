echo "Starting PoisosnSinx Regression test"
echo ${PWD}
cp ../PoissonSinx/options.input poisson_sinx_options.input
../PoissonSinx/poisson_sinx_uniform_driver poisson_sinx_options.input > disco4est.out
RES=$(cat disco4est.out | grep -c "027 FNRMSQR 0.000000000000000000000378767161")
rm poisson_sinx_options.input
if [ $RES -ne "1" ]; then
    echo "Test Fail"
    echo "****************** CDS REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** CDS REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** CDS REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** CDS REGRESSION TEST OUTPUT BEGIN *****************"
    cat disco4est.out
    echo "****************** CDS REGRESSION TEST OUTPUT END *****************"
    echo "****************** CDS REGRESSION TEST OUTPUT END *****************"
    echo "****************** CDS REGRESSION TEST OUTPUT END *****************"
    echo "****************** CDS REGRESSION TEST OUTPUT END *****************"
    rm disco4est.out
    exit 1
else
    echo "Poisson Sinx Regression Test Success"
    rm disco4est.out
fi



