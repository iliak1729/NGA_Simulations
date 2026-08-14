#!/bin/bash



# Run 64 Cases
echo "Running with input64 ..."
time mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i inputH -v 2
echo "Complete"

echo "Copying DATA to OLD file"

# cp -r ./vtk ./OLD/N64-64-96-3D-NOSHIFT/
# cp -r ./monitor ./OLD/N64-64-96-3D-NOSHIFT/

echo "COPY COMPLETE"

echo "All runs complete."