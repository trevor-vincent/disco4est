echo "Starting Stamm Regression test"
echo ${PWD}
sed 's|num_of_amr_steps = .*|num_of_amr_steps = 9|g' ../Stamm/options.input > stamm_options.input
../Stamm/stamm_driver stamm_options.input > disco4est.out
RES0=$(cat disco4est.out | grep -c "0.000000001330124")
RES1=$(cat disco4est.out | grep -c "0.00000001711726")
RES2=$(cat disco4est.out | grep -c "0.000001303205")
RES3=$(cat disco4est.out | grep -c "0.0000007629001")
RES=$(echo $RES0 + $RES1 + $RES2 + $RES3 | bc) 
rm stamm_options.input
if [ $RES -ne "4" ]; then
    echo "Test Fail"
    echo "****************** STAMM REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** STAMM REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** STAMM REGRESSION TEST OUTPUT BEGIN *****************"
    echo "****************** STAMM REGRESSION TEST OUTPUT BEGIN *****************"
    cat disco4est.out
    echo "****************** STAMM REGRESSION TEST OUTPUT END *****************"
    echo "****************** STAMM REGRESSION TEST OUTPUT END *****************"
    echo "****************** STAMM REGRESSION TEST OUTPUT END *****************"
    echo "****************** STAMM REGRESSION TEST OUTPUT END *****************"
    rm disco4est.out
    exit 1
else
    echo "Test Success"
    rm disco4est.out
fi



