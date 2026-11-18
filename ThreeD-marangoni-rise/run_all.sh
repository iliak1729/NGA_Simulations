#!/bin/bash



# Run 64 Cases
echo "Running with input64 ..."
rm -rf ./vtk
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input64_P_NoMatching -v 2
echo "Complete"

echo "Copying DATA to OLD file"
cp -r ./monitor ./NewProperties/P_NoMatching_Jibben_96/
cp -r ./vtk ./NewProperties/P_NoMatching_Jibben_96/

echo "COPY COMPLETE"

echo "All runs complete."