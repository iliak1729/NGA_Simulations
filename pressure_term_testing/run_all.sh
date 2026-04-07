#!/bin/bash



# Shift Peskin
echo "Running with inputShiftPeskin120..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputShiftPeskin120 -v 0 
echo "Complete"

echo "Running with inputShiftPeskin1200..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputShiftPeskin1200 -v 0 
echo "Complete"

echo "Running with inputShiftPeskin12000..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputShiftPeskin12000 -v 0 
echo "Complete"

# Peskin
echo "Running with inputPeskin120..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputPeskin120 -v 0 
echo "Complete"

echo "Running with inputPeskin1200..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputPeskin1200 -v 0 
echo "Complete"

echo "Running with inputPeskin12000..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputPeskin12000 -v 0 
echo "Complete"

# Seric
echo "Running with inputSeric120..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputSeric120 -v 0 
echo "Complete"

echo "Running with inputSeric1200..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputSeric1200 -v 0 
echo "Complete"

echo "Running with inputSeric12000..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputSeric12000 -v 0 
echo "Complete"


echo " Translation Cases Complete"
echo " Running Laplace Cases"


echo "Running with inputShiftPeskin120_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputShiftPeskin120_Laplace -v 0 
echo "Complete"

echo "Running with inputShiftPeskin1200_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputShiftPeskin1200_Laplace -v 0 
echo "Complete"

echo "Running with inputShiftPeskin12000_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputShiftPeskin12000_Laplace -v 0 
echo "Complete"


echo "Running with inputPeskin120_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputPeskin120_Laplace -v 0 
echo "Complete"

echo "Running with inputPeskin1200_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputPeskin1200_Laplace -v 0 
echo "Complete"

echo "Running with inputPeskin12000_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputPeskin12000_Laplace -v 0 
echo "Complete"


echo "Running with inputSeric120_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputSeric120_Laplace -v 0 
echo "Complete"

echo "Running with inputSeric1200_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputSeric1200_Laplace -v 0 
echo "Complete"

echo "Running with inputSeric12000_Laplace..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i inputSeric12000_Laplace -v 0 
echo "Complete"

echo "All runs complete."