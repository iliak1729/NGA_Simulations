#!/bin/bash



# Run 64 Cases
echo "Running with input64 ..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input64 -v 2
echo "Complete"

echo "Copying DATA to OLD file"
cp -r ./monitor ./FaceParaboloid/P_Matching_Jibben/
cp -r ./vtk ./FaceParaboloid/P_Matching_Jibben/

echo "COPY COMPLETE"

echo "All runs complete."