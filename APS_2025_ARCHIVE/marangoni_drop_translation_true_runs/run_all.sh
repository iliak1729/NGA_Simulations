#!/bin/bash

# Run CSS PLIC
echo "Running with input64..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input64 -v 0 
echo "Complete"

echo "Running with input128..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input128 -v 0 
echo "Complete"

echo "Running with input256..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input256 -v 0 
echo "Complete"

# Run CSS JIbben
echo "Running with input64P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input64P -v 0 
echo "Complete"

echo "Running with input128P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input128P -v 0 
echo "Complete"

echo "Running with input256P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input256P -v 0 
echo "Complete"


# Run 512
echo "Running with input512..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input512 -v 0 
echo "Complete"

# Run 512 Jibben
echo "Running with input512P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input512P -v 0 
echo "Complete"
