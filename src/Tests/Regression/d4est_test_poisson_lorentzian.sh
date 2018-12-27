echo "Starting Lorentzian Regression test"
echo ${PWD}
cp ../Poisson/options_lorentzian_amr.input lorentzian_options.input
../Poisson/poisson_lorentzian_amr_multigrid_anares_driver lorentzian_options.input > disco4est.out
sed -i 's|num_of_amr_steps = .*|num_of_amr_steps = 1|g'  lorentzian_options.input
RES=$(cat disco4est.out | grep -c "104 832 832 2706.02899845")
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



