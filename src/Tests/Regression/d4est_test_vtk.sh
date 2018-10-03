echo "Starting VTK Regression test"
echo ${PWD}
./d4est_test_vtk
TEST=$(cmp --silent $old $new || echo "files are different")
RES0=$(echo $TEST | grep -c different)

if [ $RES -ne "0" ]; then
    echo "d4est test vtk diff failed"
else
    echo "d4est test vtk diff passed"
fi



