# Run CSF PLIC
echo "Running with inputCSF8..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF8 -v 0 
echo "Complete"

echo "Running with inputCSF16..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF16 -v 0 
echo "Complete"

echo "Running with inputCSF32..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF32 -v 0 
echo "Complete"

echo "Running with inputCSF64..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF64 -v 0 
echo "Complete"

echo "Running with inputCSF128..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF128 -v 0 
echo "Complete"

# Run CSS PLIC
echo "Running with inputCSS8..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS8 -v 0 
echo "Complete"

echo "Running with inputCSS16..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS16 -v 0 
echo "Complete"

echo "Running with inputCSS32..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS32 -v 0 
echo "Complete"

echo "Running with inputCSS64..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS64 -v 0 
echo "Complete"

echo "Running with inputCSS128..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS128 -v 0 
echo "Complete"







