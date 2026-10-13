#!/usr/bin/env bash

# Standalone build and run test script for mizuRoute on NCAR Derecho (Intel compiler)
# Addresses ESCOMP/mizuRoute issue #647

set -e

echo "=== mizuRoute Standalone Build & Run Tester (Derecho/Intel) ==="

# 1. Environment setup
echo "Setting up module environment on Derecho..."
module purge || true
module load cmake intel cray-mpich netcdf-mpi ncarcompilers

# 2. Check externals / submodules for mizuRoute (handling both standalone and CTSM component layouts)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="$( cd "${SCRIPT_DIR}/.." && pwd )"
MIZUROUTE_ROOT="$( cd "${BUILD_DIR}/../.." && pwd )"

if [ ! -d "${MIZUROUTE_ROOT}/libraries/parallelio" ]; then
    if [ -d "${MIZUROUTE_ROOT}/../../libraries/parallelio" ]; then
        echo "Linking parallelio from parent CTSM checkout..."
        mkdir -p "${MIZUROUTE_ROOT}/libraries"
        ln -snf "${MIZUROUTE_ROOT}/../../libraries/parallelio" "${MIZUROUTE_ROOT}/libraries/parallelio"
    fi
fi

if [ ! -d "${MIZUROUTE_ROOT}/externals/toml-f" ]; then
    if [ -d "${MIZUROUTE_ROOT}/../../externals/toml-f" ]; then
        echo "Linking toml-f from parent CTSM checkout..."
        mkdir -p "${MIZUROUTE_ROOT}/externals"
        ln -snf "${MIZUROUTE_ROOT}/../../externals/toml-f" "${MIZUROUTE_ROOT}/externals/toml-f"
    fi
fi

# 3. Build standalone binary
echo "Building standalone mizuRoute binary in ${BUILD_DIR}..."
cd "${BUILD_DIR}"
export BLDDIR="${BUILD_DIR}/../"

rm -rf "${MIZUROUTE_ROOT}/externals/toml-f/_build" "${MIZUROUTE_ROOT}/externals/toml-f/_install"

gmake clean FC=intel FC_EXE=mpif90 F_MASTER="${BLDDIR}" || true
gmake FC=intel FC_EXE=mpif90 F_MASTER="${BLDDIR}" NCDF_PATH="${NETCDF}" MODE=fast EXE=route_runoff

# Verify standalone executable binary was created in route/bin/ by gmake install
EXE_BIN="${MIZUROUTE_ROOT}/route/bin/route_runoff"
if [ ! -f "${EXE_BIN}" ]; then
    echo "ERROR: Standalone executable 'route_runoff' was not created at ${EXE_BIN}."
    exit 1
fi

echo "SUCCESS: Standalone executable 'route_runoff' built successfully at ${EXE_BIN}."
echo "=== Standalone Build Test Passed ==="
