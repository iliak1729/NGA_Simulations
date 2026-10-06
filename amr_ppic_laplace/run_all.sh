#!/bin/bash

CASES=(
    "input_120CSF"
)

# Core counts to test
CORES=(
    8
    12
    16
)

for INPUT in "${CASES[@]}"; do

    # Remove "input_" from the name for the output directory
    CASE="${INPUT#input_}"

    for NCORES in "${CORES[@]}"; do

        DIR="./OLD_LapSweep_PPIC_ChanginCores/amr_LA${CASE}/${NCORES}Cores"

        echo "========================================"
        echo "Running case: $CASE"
        echo "Input file:   $INPUT"
        echo "Cores:        $NCORES"
        echo "Output dir:   $DIR"
        echo "========================================"

        echo "Removing old Viz Data"
        rm -rf ./amrviz/

        echo "Running $INPUT on $NCORES cores ..."
        time mpiexec -n "$NCORES" ./nga2.dp.gnu.opt.mpi.exe -i "$INPUT" -v 0

        echo "Saving monitor data"
        mkdir -p "$DIR"
        cp -r ./monitor "$DIR/"
        echo "Saving amrviz data"
        cp -r ./amrviz "$DIR/"

        echo "Case $CASE on $NCORES cores complete"
        echo

    done

done

echo "All simulations complete."