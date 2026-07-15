#!/bin/bash


# Run 128 Cases
# echo "Running with input2D ..."
# time mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input2D -v 0
# echo "Complete"

# echo "Copying Data"
# cp ./monitor/Error_Values ./monitor/Error_Values_2D_ProjectedMetric_MeanCurvatureSweep
# echo "COPY COMPLETE"

echo "Running with input3D ..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i input3D -v 0
echo "Complete"

echo "Copying Data"
cp ./monitor/Error_Values ./monitor/Error_Values_3D_ProjectedMetric_MeanCurvatureSweep
echo "COPY COMPLETE"

echo "All runs complete."