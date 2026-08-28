#!/bin/bash



# Run 64 Cases
echo "Removing old Viz Data"
rm -rf ./amrviz/
echo "Running with input ..."
time mpiexec -n 24 ./nga2.dp.gnu.opt.mpi.exe -i input -v 2
echo "Complete"

echo "Copying DATA to OLD file"

# cp -r ./vtk ./OLD/N64-64-96-3D-NOSHIFT/
# cp -r ./monitor ./OLD/N64-64-96-3D-NOSHIFT/

echo "COPY COMPLETE"

echo "All runs complete."