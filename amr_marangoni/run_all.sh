
CASES=(
    "input_Level5P"
)

for INPUT in "${CASES[@]}"; do

    # Remove "input_" from the name for the output directory
    CASE="${INPUT#input_}"
    DIR="./APS2026FINAL/amr_${CASE}"

    echo "========================================"
    echo "Running case: $CASE"
    echo "Input file:   $INPUT"
    echo "Output dir:   $DIR"
    echo "========================================"

    echo "Removing old Viz Data"
    rm -rf ./amrviz/

    echo "Running with $INPUT ..."
    time mpiexec -n 8 ./nga2.dp.gnu.opt.mpi.exe -i "$INPUT" -v 0

    echo "Saving monitor data"
    mkdir -p "$DIR"
    cp -r ./monitor "$DIR/"

    echo "Case $CASE complete"
    echo

done

echo "All simulations complete."