#!/usr/bin/env bash

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Usage: ./Allrun.sh [initial-cell-count] [number-of-refinements]
GRID_CELLS="${1:-100}"
REFINEMENTS="${2:-19}"

[[ "$GRID_CELLS" =~ ^[1-9][0-9]*$ ]] || {
    echo "Initial cell count must be a positive integer: $GRID_CELLS" >&2
    exit 2
}
[[ "$REFINEMENTS" =~ ^[0-9]+$ ]] || {
    echo "Number of refinements must be a non-negative integer: $REFINEMENTS" >&2
    exit 2
}

set +e +u
source "$HOME/SoftwareBuilds/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
OPENFOAM_SOURCE_STATUS=$?
if (( OPENFOAM_SOURCE_STATUS != 0 )); then
    echo "Failed to source the OpenFOAM environment (status $OPENFOAM_SOURCE_STATUS)." >&2
    exit "$OPENFOAM_SOURCE_STATUS"
fi

cd "$CASE_DIR"
set -euo pipefail

INITIAL_FIELDS_DIR="$(mktemp -d)"
cp -a 0/. "$INITIAL_FIELDS_DIR/"
trap 'rm -rf -- "$INITIAL_FIELDS_DIR"' EXIT

backup_existing_results() {
    local results_dir="$1"
    local suffix

    if [[ -d "$results_dir" ]]; then
        suffix="$(date +%Y%m%d_%H%M%S)"
        mv -- "$results_dir" "${results_dir}_preexisting_${suffix}"
    fi
}

clean_generated_case() {
    local entry name

    rm -rf -- constant/polyMesh

    # Keep the initial fields in 0; all solver-produced time directories are regenerated.
    for entry in ./*; do
        [[ -d "$entry" ]] || continue
        name="${entry#./}"
        if [[ "$name" != 0 && "$name" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
            rm -rf -- "$entry"
        fi
    done
}

set_mesh_cell_count() {
    local cells="$1"

    sed -E -i \
        "s/^([[:space:]]*)\\([0-9]+[[:space:]]+1[[:space:]]+1\\)([[:space:]]*\\/\\/.*)?$/\\1(${cells} 1 1)\\2/" \
        system/blockMeshDict
}

reset_initial_fields() {
    rm -rf -- 0
    mkdir -- 0
    cp -a "$INITIAL_FIELDS_DIR"/. 0/

    # Use uniform templates so setFields can create fields for any refined mesh size.
    foamDictionary 0/p -entry internalField -set "uniform 100000"
    foamDictionary 0/T -entry internalField -set "uniform 300"
    foamDictionary 0/U -entry internalField -set "uniform (0 0 0)"
}

run_verification() {
    local cells="$1"
    local result_dir="verificationResults_${cells}"

    echo
    echo "============================================================"
    echo "Running normal-shock verification with ${cells} cells"
    echo "============================================================"

    set_mesh_cell_count "$cells"
    clean_generated_case
    reset_initial_fields
    backup_existing_results verificationResults

    blockMesh | tee log.blockMesh
    checkMesh | tee log.checkMesh
    setFields | tee log.setFields
    rhoCentralFoam | tee log.rhoCentralFoam

    touch normalShockVerification.foam
    postProcess -func writeCellCentres -latestTime

    CASE_DIR="$CASE_DIR" OUTPUT_DIR="$CASE_DIR/verificationResults" \
        julia scripts/verifyNormalShock.jl | tee log.verifyNormalShock

    if [[ -e "$result_dir" ]]; then
        result_dir="${result_dir}_$(date +%Y%m%d_%H%M%S)"
    fi
    mv -- verificationResults "$result_dir"
    cp -- log.blockMesh log.checkMesh log.setFields log.rhoCentralFoam \
        log.verifyNormalShock "$result_dir/"

    echo "Saved verification results to $result_dir"
}

for ((run = 0; run <= REFINEMENTS; run++)); do
    run_verification "$GRID_CELLS"

    # Round to the nearest integer so every refinement grows by approximately 20%.
    if (( run < REFINEMENTS )); then
        GRID_CELLS=$(( (GRID_CELLS * 6 + 5) / 5 ))
    fi
done

julia scripts/gridRefinementStudy.jl >> log.gridRefinementStudy