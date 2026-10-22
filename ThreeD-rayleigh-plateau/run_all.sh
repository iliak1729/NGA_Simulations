
CASES=(
    "input_PCST_P_Matching"
)
make -j 10
for INPUT in "${CASES[@]}"; do

    # Remove "input_" from the name for the output directory
    CASE="${INPUT#input_}"
    DIR="./FACE_PARABOLOID/Debugging_${CASE}"

    echo "========================================"
    echo "Running case: $CASE"
    echo "Input file:   $INPUT"
    echo "Output dir:   $DIR"
    echo "========================================"

    echo "Removing old Viz Data"
    rm -rf ./vtk/

    echo "Running with $INPUT ..."
    time mpiexec -n 16 ./nga.dp.gnu.opt.mpi.exe -i "$INPUT" -v 0
    
    echo "Saving monitor data"
    mkdir -p "$DIR"
    cp -r ./monitor "$DIR/"
    cp -r ./vtk "$DIR/"

    echo "Case $CASE complete"
    echo

done

echo "All simulations complete."