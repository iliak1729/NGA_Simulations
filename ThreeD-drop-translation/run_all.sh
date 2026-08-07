#!/bin/bash



# Run 2D Cases
# echo "Running with input2D_N64_L120_PUST_LVIRA ..."
# time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i input2D_N64_L120_PUST_LVIRA -v 0
# echo "Complete"

# echo "Copying DATA to OLD file"

# mkdir -p ./OLD/input2D_N64_L120_PUST_LVIRA/
# cp -r ./monitor ./OLD/input2D_N64_L120_PUST_LVIRA/

# echo "COPY COMPLETE"




echo "Running with input2D_N64_L120_PUST_Jibben ..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i input2D_N64_L120_PUST_Jibben -v 0
echo "Complete"

echo "Copying DATA to OLD file"

mkdir -p ./OLD/input2D_N64_L120_PUST_Jibben/
cp -r ./monitor ./OLD/input2D_N64_L120_PUST_Jibben/

echo "COPY COMPLETE"




echo "Running with input2D_N64_L1200_PUST_LVIRA ..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i input2D_N64_L1200_PUST_LVIRA -v 0
echo "Complete"

echo "Copying DATA to OLD file"
mkdir -p ./OLD/input2D_N64_L1200_PUST_LVIRA/
cp -r ./monitor ./OLD/input2D_N64_L1200_PUST_LVIRA/

echo "COPY COMPLETE"




echo "Running with input2D_N64_L1200_PUST_Jibben ..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i input2D_N64_L1200_PUST_Jibben -v 0
echo "Complete"

echo "Copying DATA to OLD file"

mkdir -p ./OLD/input2D_N64_L1200_PUST_Jibben/
cp -r ./monitor ./OLD/input2D_N64_L1200_PUST_Jibben/

echo "COPY COMPLETE"




echo "Running with input2D_N64_L12000_PUST_LVIRA ..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i input2D_N64_L12000_PUST_LVIRA -v 0
echo "Complete"

echo "Copying DATA to OLD file"
mkdir -p ./OLD/input2D_N64_L12000_PUST_LVIRA/
cp -r ./monitor ./OLD/input2D_N64_L12000_PUST_LVIRA/

echo "COPY COMPLETE"




echo "Running with input2D_N64_L12000_PUST_Jibben ..."
time mpiexec -n 9 ./nga.dp.gnu.opt.mpi.exe -i input2D_N64_L12000_PUST_Jibben -v 0
echo "Complete"

echo "Copying DATA to OLD file"
mkdir -p ./OLD/input2D_N64_L12000_PUST_Jibben/
cp -r ./monitor ./OLD/input2D_N64_L12000_PUST_Jibben/

echo "COPY COMPLETE"

echo "All runs complete."