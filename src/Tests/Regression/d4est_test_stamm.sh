echo "Starting Stamm Regression test"
echo ${PWD}
sed 's|num_of_amr_steps = .*|num_of_amr_steps = 9|g' ../Stamm/options.input > stamm_options.input
# RES=$(../Stamm/stamm_driver stamm_options.input | grep -c "0.0000001165842921684272017 0.0000007074379709713716485 0.0000146083448435542040227 0.0001316349774934072013998")
../Stamm/stamm_driver stamm_options.input > disco4est.out
# RES=$(cat disco4est.out | grep -c "0.0000000013301248793227330 0.0000000171172621537055147 0.0000013032058553497301960 0.0000007629001363926490040")
RES0=$(cat disco4est.out | grep -c "0.000000001330124879")
RES1=$(cat disco4est.out | grep -c "0.000000017117262153")
RES2=$(cat disco4est.out | grep -c "0.0000013032058553")
RES3=$(cat disco4est.out | grep -c "0.000000762900136392")
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



