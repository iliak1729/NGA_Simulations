CASES=(
    # ========================================
    # D8
    # ========================================
    Basak_input_CSF_D8_SingleFluid
    Basak_input_PCST_D8_NoP_NoMatching_SingleFluid
    Basak_input_PCST_D8_NoP_Matching_SingleFluid
    Basak_input_PCST_D8_P_NoMatching_SingleFluid
    Basak_input_PCST_D8_P_Matching_SingleFluid

    # Basak_input_CSF_D8_DiffFluid
    # Basak_input_PCST_D8_NoP_NoMatching_DiffFluid
    # Basak_input_PCST_D8_NoP_Matching_DiffFluid
    # Basak_input_PCST_D8_P_NoMatching_DiffFluid
    # Basak_input_PCST_D8_P_Matching_DiffFluid

    # ========================================
    # D16
    # ========================================
    Basak_input_CSF_D16_SingleFluid
    Basak_input_PCST_D16_NoP_NoMatching_SingleFluid
    Basak_input_PCST_D16_NoP_Matching_SingleFluid
    Basak_input_PCST_D16_P_NoMatching_SingleFluid
    Basak_input_PCST_D16_P_Matching_SingleFluid

    # Basak_input_CSF_D16_DiffFluid
    # Basak_input_PCST_D16_NoP_NoMatching_DiffFluid
    # Basak_input_PCST_D16_NoP_Matching_DiffFluid
    # Basak_input_PCST_D16_P_NoMatching_DiffFluid
    # Basak_input_PCST_D16_P_Matching_DiffFluid

    # ========================================
    # D32
    # ========================================
    Basak_input_CSF_D32_SingleFluid
    Basak_input_PCST_D32_NoP_NoMatching_SingleFluid
    Basak_input_PCST_D32_NoP_Matching_SingleFluid
    Basak_input_PCST_D32_P_NoMatching_SingleFluid
    Basak_input_PCST_D32_P_Matching_SingleFluid

    # Basak_input_CSF_D32_DiffFluid
    # Basak_input_PCST_D32_NoP_NoMatching_DiffFluid
    # Basak_input_PCST_D32_NoP_Matching_DiffFluid
    # Basak_input_PCST_D32_P_NoMatching_DiffFluid
    # Basak_input_PCST_D32_P_Matching_DiffFluid

    # ========================================
    # D64
    # ========================================
    # Basak_input_CSF_D64_SingleFluid
    # Basak_input_PCST_D64_NoP_NoMatching_SingleFluid
    # Basak_input_PCST_D64_NoP_Matching_SingleFluid
    # Basak_input_PCST_D64_P_NoMatching_SingleFluid
    # Basak_input_PCST_D64_P_Matching_SingleFluid

    # Basak_input_CSF_D64_DiffFluid
    # Basak_input_PCST_D64_NoP_NoMatching_DiffFluid
    # Basak_input_PCST_D64_NoP_Matching_DiffFluid
    # Basak_input_PCST_D64_P_NoMatching_DiffFluid
    # Basak_input_PCST_D64_P_Matching_DiffFluid
)
make -j 10
for INPUT in "${CASES[@]}"; do

    # Remove "input_" from the name for the output directory
    CASE="${INPUT#Basak_input_}"
    DIR="./Basak_Jibben/StaticDrop_${CASE}"

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