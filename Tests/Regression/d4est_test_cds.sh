echo "Starting ConstantDensityStar Regression test"
echo ${PWD}
sed 's|num_of_amr_steps = .*|num_of_amr_steps = 1|g' ../ConstantDensityStar/options.input > options.input
sed -i 's|ksp_type = .*|ksp_type = cg|g' options.input
sed -i 's|ksp_do_not_use_preconditioner = 0|ksp_do_not_use_preconditioner = 1|g' options.input
../ConstantDensityStar/constant_density_star_driver options.input 2&>1 disco4est.out
RES=$(cat disco4est.out | grep -c 0.0000318356651187705076136)
rm options.input
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
    echo "CDS Regression Test Success"
    rm disco4est.out
fi



