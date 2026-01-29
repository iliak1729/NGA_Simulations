#!/bin/bash



# CSS PLIC
echo "Running with inputCSS120..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputCSS120 -v 0 
echo "Complete"

echo "Running with inputCSS1200..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputCSS1200 -v 0 
echo "Complete"

echo "Running with inputCSS12000..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputCSS12000 -v 0 
echo "Complete"

echo "All runs complete."