
CASES=(
    "input_PCST_P_NoMatching_BothRatio"
)
make -j 10
for INPUT in "${CASES[@]}"; do

    # Remove "input_" from the name for the output directory
    CASE="${INPUT#input_}"
    # DIR="./Debugging_PostOvernight/JibbenSq_ShortTime_${CASE}"
    DIR="./APS2026FINAL_DEBUG/Domain64_${CASE}"

    echo "========================================"
    echo "Running case: $CASE"
    echo "Input file:   $INPUT"
    echo "Output dir:   $DIR"
    echo "========================================"

    echo "Removing old Viz Data"
    rm -rf ./vtk/

    echo "Running with $INPUT ..."
    time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i "$INPUT" -v 2
    
    echo "Saving monitor data"
    mkdir -p "$DIR"
    cp -r ./monitor "$DIR/"
    cp -r ./vtk "$DIR/"

    echo "Case $CASE complete"
    echo

done

echo "All simulations complete."