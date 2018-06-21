echo "Starting Lorentzian Regression test"
echo ${PWD}
cp ../Poisson/options_lorentzian_mg_test.input lorentzian_options.input
../Poisson/poisson_lorentzian_multigrid_driver lorentzian_options.input > disco4est.out
RES=$(cat disco4est.out | grep -c "8 216 216 2.3466279072")
rm lorentzian_options.input
if [ $RES -ne "1" ]; then
    echo "Test Fail"
    echo "****************** LORENTZIAN REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** LORENTZIAN REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** LORENTZIAN REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** LORENTZIAN REGRESSION TEST OUTPUT BEGIN *****************"
    cat disco4est.out
    echo "****************** LORENTZIAN REGRESSION TEST OUTPUT END *****************"
    echo "****************** LORENTZIAN REGRESSION TEST OUTPUT END *****************"
    echo "****************** LORENTZIAN REGRESSION TEST OUTPUT END *****************"
    echo "****************** LORENTZIAN REGRESSION TEST OUTPUT END *****************"
    rm disco4est.out
    exit 1
else
    echo "Lorentzian Regression Test Success"
    rm disco4est.out
fi



