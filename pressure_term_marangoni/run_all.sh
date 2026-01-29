#!/bin/bash

# Run CSS PLIC
echo "Running with input32..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input32 -v 0 
echo "Complete"

echo "Running with input64..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input64 -v 0 
echo "Complete"

echo "Running with input128..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input128 -v 0 
echo "Complete"

echo "All runs complete."