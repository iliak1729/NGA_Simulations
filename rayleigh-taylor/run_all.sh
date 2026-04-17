#!/bin/bash

# Run 64 Cases
# echo "Running with0 inputSeric64 ..."
# time mpiexec -n 6 ./nga.dp.gnu.opt.mpi.exe -i inputSeric64 -v 0 
# echo "Complete"

# echo "Running with inputPeskin64 ..."
# time mpiexec -n 6 ./nga.dp.gnu.opt.mpi.exe -i inputPeskin64 -v 0 
# echo "Complete"

# echo "Running with inputShiftPeskin64 ..."
# time mpiexec -n 6 ./nga.dp.gnu.opt.mpi.exe -i inputShiftPeskin64 -v 0 
# echo "Complete"

# Run 128 Cases
echo "Running with0 inputSeric128 ..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputSeric128 -v 0 
echo "Complete"

echo "Running with inputPeskin128 ..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputPeskin128 -v 0 
echo "Complete"

echo "Running with inputShiftPeskin128 ..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputShiftPeskin128 -v 0 
echo "Complete"

echo "All runs complete."