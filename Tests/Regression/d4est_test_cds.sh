echo "Starting ConstantDensityStar Regression test"
echo ${PWD}
sed 's|num_of_amr_steps = .*|num_of_amr_steps = 1|g' ../ConstantDensityStar/options.input > options.input
sed -i 's|ksp_type = .*|ksp_type = cg|g' options.input
sed -i 's|ksp_do_not_use_preconditioner = 0|ksp_do_not_use_preconditioner = 1|g' options.input
RES=$(../ConstantDensityStar/constant_density_star_driver options.input | grep -c 0.0000318356651187705076136)
rm options.input
if [ $RES -ne "1" ]; then
    echo "Test Fail"
    exit 1
else
    echo "Test Success"
fi



