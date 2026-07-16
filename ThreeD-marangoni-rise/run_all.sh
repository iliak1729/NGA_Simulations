#!/bin/bash


# Run 128 Cases
echo "Running with input128 ..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input128 -v 2
echo "Complete"

echo "Copying DATA to OLD file"

cp -r ./vtk ./OLD/N128-128-192-3D/
cp -r ./monitor ./OLD/N128-128-192-3D/

echo "COPY COMPLETE"

echo "All runs complete."