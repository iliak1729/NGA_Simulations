#!/bin/bash


# CSF PLIC
echo "Running with inputCSF120..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSF120 -v 0 
echo "Complete"

echo "Running with inputCSF1200..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSF1200 -v 0 
echo "Complete"

echo "Running with inputCSF12000..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSF12000 -v 0 
echo "Complete"

# CSS PLIC
echo "Running with inputCSS120..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSS120 -v 0 
echo "Complete"

echo "Running with inputCSS1200..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSS1200 -v 0 
echo "Complete"

echo "Running with inputCSS12000..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSS12000 -v 0 
echo "Complete"

# CSF Jibben
echo "Running with inputCSF120P..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSF120P -v 0 
echo "Complete"

echo "Running with inputCSF1200P..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSF1200P -v 0 
echo "Complete"

echo "Running with inputCSF12000P..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSF12000P -v 0 
echo "Complete"

# CSS Jibben
echo "Running with inputCSS120P..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSS120P -v 0 
echo "Complete"

echo "Running with inputCSS1200P..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSS1200P -v 0 
echo "Complete"

echo "Running with inputCSS12000P..."
time mpiexec -n 24 ./nga.dp.gnu.opt.mpi.exe -i inputCSS12000P -v 0 
echo "Complete"




