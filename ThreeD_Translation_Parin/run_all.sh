CASES=(
    # ========================================
    # D8
    # ========================================
    "input_La120_We04_D8"
    "input_La1200_We04_D8"
    "input_La12000_We04_D8"
)
make -j 10
for INPUT in "${CASES[@]}"; do

    # Remove "input_" from the name for the output directory
    CASE="${INPUT#Basak_input_}"
    DIR="./Final2026APS_Translation/JibbenPU_${CASE}"
    # DIR="./Debugging/StaticDrop_${CASE}"

    echo "========================================"
    echo "Running case: $CASE"
    echo "Input file:   $INPUT"
    echo "Output dir:   $DIR"
    echo "========================================"

    echo "Removing old Viz Data"
    rm -rf ./vtk/

    echo "Running with $INPUT ..."
    echo "Simulation start time: $(date '+%Y-%m-%d %H:%M:%S')"
    time mpiexec -n 8 ./nga.dp.gnu.opt.mpi.exe -i "$INPUT" -v 0
    
    echo "Saving monitor data"
    mkdir -p "$DIR"
    cp -r ./monitor "$DIR/"
    cp -r ./vtk "$DIR/"

    echo "Case $CASE complete"
    echo

done

echo "All simulations complete."