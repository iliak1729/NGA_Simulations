#!/bin/bash



# CSS PLIC
echo "Running with inputCSS16..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputCSS16 -v 0 
echo "Complete"

echo "Running with inputCSS16_Poisson..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputCSS16_Poisson -v 0 
echo "Complete"

echo "Running with inputCSS32..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputCSS32 -v 0 
echo "Complete"

echo "Running with inputCSS32_Poisson..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputCSS32_Poisson -v 0 
echo "Complete"

echo "Running with inputCSS64..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputCSS64 -v 0 
echo "Complete"

echo "Running with inputCSS64_Poisson..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputCSS64_Poisson -v 0 
echo "Complete"





echo "All runs complete."