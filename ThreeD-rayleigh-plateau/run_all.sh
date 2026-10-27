
CASES=(
    # CSF Cases
    "Basak_input_CSF_LVIRA"
    "Basak_input_CSF_PLICNET"

    # PCST - No Pressure, Matching
    "Basak_input_PCST_NoP_M_LVIRA"
    "Basak_input_PCST_NoP_M_PLICNET"

    # PCST - No Pressure, No Matching
    "Basak_input_PCST_NoP_NoM_LVIRA"
    "Basak_input_PCST_NoP_NoM_PLICNET"

    # PCST - Pressure, Matching
    "Basak_input_PCST_P_M_LVIRA"
    "Basak_input_PCST_P_M_PLICNET"

    # PCST - Pressure, No Matching
    "Basak_input_PCST_P_NoM_LVIRA"
    "Basak_input_PCST_P_NoM_PLICNET"

    # Jibben Cases
    "Basak_input_CSF_Jibben"
    "Basak_input_PCST_NoP_M_Jibben"
    "Basak_input_PCST_NoP_NoM_Jibben"
    "Basak_input_PCST_P_M_Jibben"
    "Basak_input_PCST_P_NoM_Jibben"

    # "input_PCST_P_Matching"
    # "input_PCST_NoP_Matching"
    # "input_PCST_NoP_NoMatching"
    # "input_PCST_P_NoMatching"
    # "input_PCST_P_NoMatching_LVIRA"
    # "input_PCST_P_NoMatching_PLICNET"
    # "input_PCST_P_NoMatching_Jibben"
    # "input_CSF_PLICNET"
    # "input_CSF_LVIRA"
    # "input_CSF_Jibben"
    # "input_PCST_P_Matching_Jibben"
    # "input_PCST_P_Matching_LVIRA"
    # "input_PCST_P_Matching_PLICNET"
    # "input_CSF"
    # "input_PCST_NoP_Matching_Laplacian"
    # "input_PCST_NoP_NoMatching_Laplacian"
    # "input_PCST_P_NoMatching_Laplacian"
    # "input_PCST_P_Matching_Laplacian"
)
make -j 10
for INPUT in "${CASES[@]}"; do

    # Remove "input_" from the name for the output directory
    CASE="${INPUT#input_}"
    DIR="./Basak/First_${CASE}"

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