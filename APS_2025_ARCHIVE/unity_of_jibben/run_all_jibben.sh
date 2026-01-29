# Run CSS Jibben
echo "Running with inputCSS8P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS8P -v 0 
echo "Complete"

echo "Running with inputCSS16P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS16P -v 0 
echo "Complete"

echo "Running with inputCSS32P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS32P -v 0 
echo "Complete"

echo "Running with inputCSS64P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS64P -v 0 
echo "Complete"

echo "Running with inputCSS128P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSS128P -v 0 
echo "Complete"

# Run CSF Jibben
echo "Running with inputCSF8P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF8P -v 0 
echo "Complete"

echo "Running with inputCSF16P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF16P -v 0 
echo "Complete"

echo "Running with inputCSF32P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF32P -v 0 
echo "Complete"

echo "Running with inputCSF64P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF64P -v 0 
echo "Complete"

echo "Running with inputCSF128P..."
time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i inputCSF128P -v 0 
echo "Complete"