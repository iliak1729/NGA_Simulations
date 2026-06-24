#!/bin/bash


# Run 128 Cases
echo "Running with input64 ..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input64 -v 2
echo "Complete"

echo "Copying DATA to OLD file"

cp -r ./vtk ./OLD/N64-64-96-LowCP/
cp -r ./monitor ./OLD/N64-64-96-LowCP/

echo "COPY COMPLETE"

echo "All runs complete."