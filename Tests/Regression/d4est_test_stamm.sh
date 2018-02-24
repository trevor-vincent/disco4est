#1/bin/bash

echo "Starting Stamm Regression test"
echo ${PWD}
sed -i 's|num_of_amr_steps = .*|num_of_amr_steps = 4|g' ./Stamm/options.input
RES=$(./Stamm/stamm_driver ./Stamm/options.input | grep -c "0.0000001165842921684272017 0.0000007074379709713716485 0.0000146083448435542040227 0.0001316349774934072013998")
if [ $RES -ne "1" ]; then
    echo "Test Fail"
    exit 1
else
    echo "Test Success"
fi



