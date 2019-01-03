echo "Starting Schwarz Regression test"
echo ${PWD}
../Tests/d4est_test_schwarz_cubed_sphere_new_3d > disco4est.out
# sed -i 's|num_of_amr_steps = .*|num_of_amr_steps = 1|g'  lorentzian_options.input
RES=$(cat disco4est.out | grep -c "0.15228638")
rm lorentzian_options.input
if [ $RES -ne "1" ]; then
    echo "Test Fail"
    echo "****************** SCHWARZ REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** SCHWARZ REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** SCHWARZ REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** SCHWARZ REGRESSION TEST OUTPUT BEGIN *****************"
    cat disco4est.out
    echo "****************** SCHWARZ REGRESSION TEST OUTPUT END *****************"
    echo "****************** SCHWARZ REGRESSION TEST OUTPUT END *****************"
    echo "****************** SCHWARZ REGRESSION TEST OUTPUT END *****************"
    echo "****************** SCHWARZ REGRESSION TEST OUTPUT END *****************"
    rm disco4est.out
    exit 1
else
    echo "Schwarz Regression Test Success"
    rm disco4est.out
fi



