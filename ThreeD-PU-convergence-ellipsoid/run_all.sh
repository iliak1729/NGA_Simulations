#!/bin/bash


# Run 2D
echo "Running with input2D ..."
time mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input2D -v 0
echo "Complete"

echo "Copying Data"
cp ./monitor/Error_Values ./monitor/Error_Values_2D_EllipsoidStress_SpreadSweepACOS_NewNeighborhoodTest
echo "COPY COMPLETE"
# Run 3D
echo "Running with input3D ..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input3D -v 0
echo "Complete" 

echo "Copying Data"
cp ./monitor/Error_Values ./monitor/Error_Values_3D_EllipsoidStress_SpreadSweepACOS_NewNeighborhoodTest
echo "COPY COMPLETE"
# Run 2D Jibben
echo "Running with input2D ..."
time mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input2D_Jibben -v 0
echo "Complete"

echo "Copying Data"
cp ./monitor/Error_Values ./monitor/Error_Values_2D_EllipsoidStress_SpreadSweepACOS_NewNeighborhoodTest_Jibben
echo "COPY COMPLETE"
# Run 3D Jibben
echo "Running with input3D ..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input3D_Jibben -v 0
echo "Complete" 

echo "Copying Data"
cp ./monitor/Error_Values ./monitor/Error_Values_3D_EllipsoidStress_SpreadSweepACOS_NewNeighborhoodTest_Jibben
echo "COPY COMPLETE"

echo "All runs complete."