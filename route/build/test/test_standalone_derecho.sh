#!/usr/bin/env bash

# Standalone build and run test script for mizuRoute on NCAR Derecho
# Addresses ESCOMP/mizuRoute issue #647

set -e

TARGET_COMPILER="${1:-intel}"

# 1. Environment setup
if type module >/dev/null 2>&1; then
    case "${TARGET_COMPILER}" in
        intel|Intel)
            FC_COMPILER="intel"
            FC_EXE="mpif90"
            MODULE_LIST="cmake intel cray-mpich netcdf ncarcompilers"
            ;;
        gnu|gcc|GNU)
            FC_COMPILER="gnu"
            FC_EXE="mpif90"
            MODULE_LIST="cmake gcc cray-mpich netcdf ncarcompilers"
            ;;
        *)
            echo "Unknown compiler '${TARGET_COMPILER}'. Choose from: intel, gnu."
            exit 1
            ;;
    esac
    echo "Setting up module environment on Derecho for ${FC_COMPILER} (${FC_EXE})..."
    module --force purge || true
    module load ncarenv/25.10 ${MODULE_LIST}
    NCDF_DIR="${NETCDF}"
else
    echo "ERROR: 'module' command not found. This script is intended to run on NCAR Derecho."
    exit 1
fi

# Select GNU Make: 'gmake' on BSD/macOS (Homebrew), 'make' on standard Linux
if command -v gmake >/dev/null 2>&1; then
    MAKE_CMD="gmake"
else
    MAKE_CMD="make"
fi

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

if [ ! -d "${MIZUROUTE_ROOT}/libraries/toml-f" ]; then
    if [ -d "${MIZUROUTE_ROOT}/../../libraries/toml-f" ]; then
        echo "Linking toml-f from parent CTSM checkout..."
        mkdir -p "${MIZUROUTE_ROOT}/libraries"
        ln -snf "${MIZUROUTE_ROOT}/../../libraries/toml-f" "${MIZUROUTE_ROOT}/libraries/toml-f"
    fi
fi

# 3. Build standalone binary
echo "Building standalone mizuRoute binary in ${BUILD_DIR}..."
cd "${BUILD_DIR}"
export BLDDIR="${BUILD_DIR}/../"

rm -rf "${MIZUROUTE_ROOT}/libraries/toml-f/_build" "${BUILD_DIR}/lib/tomllib"

# Clean previous build and libraries using built-in Makefile targets (clean, cleanlibs)
"${MAKE_CMD}" clean cleanlibs FC="${FC_COMPILER}" FC_EXE="${FC_EXE}" F_MASTER="${BLDDIR}" || true
"${MAKE_CMD}" FC="${FC_COMPILER}" FC_EXE="${FC_EXE}" F_MASTER="${BLDDIR}" NCDF_PATH="${NCDF_DIR}" MODE=fast EXE=route_runoff
rm -rf "${BUILD_DIR}/route_runoff.dSYM"

# Verify standalone executable binary was created in route/bin/ by gmake install
EXE_BIN="${MIZUROUTE_ROOT}/route/bin/route_runoff"
if [ ! -f "${EXE_BIN}" ]; then
    echo "ERROR: Standalone executable 'route_runoff' was not created at ${EXE_BIN}."
    exit 1
fi

echo "SUCCESS: Standalone executable 'route_runoff' built successfully at ${EXE_BIN}."

# 4. Run standalone smoke test using SAMPLE_toml
echo "Setting up standalone smoke test case from SAMPLE_toml..."
RUN_DIR="$(mktemp -d)"
trap 'rm -rf "${RUN_DIR}"' EXIT

# Copy hydraulic routing parameter namelist template from repository default (not created from scratch)
cp "${MIZUROUTE_ROOT}/route/ancillary_data/param.nml.default" "${RUN_DIR}/param.nml"

# Minimal 1-reach synthetic river topology via inline CDL:
# - seg_id=1, tosegment=0 (drains out of basin), Basin_Area=1 km^2
# - hruid=1, seg_hru_id=1 (Hydrologic Response Unit mapping catchment to river reach)
cat << 'EOF' > "${RUN_DIR}/ntopo.cdl"
netcdf ntopo {
dimensions:
    seg = 1 ;
    hru = 1 ;
variables:
    int seg_id(seg) ;
    int tosegment(seg) ;
    double Length(seg) ;
    double Slope(seg) ;
    int hruid(hru) ;
    int seg_hru_id(hru) ;
    double Basin_Area(hru) ;
data:
    seg_id = 1 ;
    tosegment = 0 ;
    Length = 1000.0 ;
    Slope = 0.001 ;
    hruid = 1 ;
    seg_hru_id = 1 ;
    Basin_Area = 1000000.0 ;
}
EOF
ncgen -b -o "${RUN_DIR}/ntopo.nc" "${RUN_DIR}/ntopo.cdl"

# Minimal 2-day synthetic runoff dataset via inline CDL:
# - hruid: Hydrologic Response Unit identifier matching topology seg_hru_id
# - qsim: Simulated runoff flux (0.001 mm/s) generated for HRU 1 across 2 time steps
cat << 'EOF' > "${RUN_DIR}/runoff.cdl"
netcdf runoff {
dimensions:
    time = 2 ;
    hru = 1 ;
variables:
    double time(time) ;
        time:units = "days since 2000-01-01 00:00:00" ;
        time:calendar = "noleap" ;
    int hruid(hru) ;
    double qsim(time, hru) ;
        qsim:units = "mm/s" ;
        qsim:_FillValue = -9999.0 ;
data:
    time = 0.0, 1.0 ;
    hruid = 1 ;
    qsim = 0.001, 0.001 ;
}
EOF
ncgen -b -o "${RUN_DIR}/runoff.nc" "${RUN_DIR}/runoff.cdl"

# Configure mizuRoute.toml from SAMPLE_toml:
# - is_remap=false: Runoff and network share 1:1 topology, no ESMF mapping weight file needed
# - seg_id=-9999: Sentinel indicating all river reaches should be routed (no reach subsetting)
# - param.nml: Parameter namelist copied from repo default above
sed \
    -e 's|"CASE_NAME"|"standalone_smoke"|' \
    -e 's|sim_start = "yyyy-mm-dd hh:mm:ss"|sim_start = "2000-01-01 00:00:00"|' \
    -e 's|sim_end = "yyyy-mm-dd hh:mm:ss"|sim_end = "2000-01-02 00:00:00"|' \
    -e 's|"INPUT_RESTART_NC"|"coldstart"|' \
    -e "s|\"ANCIL_DIR\"|\"${RUN_DIR}/\"|" \
    -e "s|\"INPUT_DIR\"|\"${RUN_DIR}/\"|" \
    -e "s|\"OUTPUT_DIR\"|\"${RUN_DIR}/\"|" \
    -e 's|"NTOPO_NC"|"ntopo.nc"|' \
    -e 's|"DIMNAME_SEG"|"seg"|' \
    -e 's|"DIMNAME_HRU"|"hru"|' \
    -e 's|"RUNOFF_NC"|"runoff.nc"|' \
    -e 's|"VARNAME_RUNOFF"|"qsim"|' \
    -e 's|"VARNAME_TIME"|"time"|' \
    -e 's|"VARNAME_RO_HRU"|"hruid"|' \
    -e 's|"DIMNAME_TIME"|"time"|' \
    -e 's|"DIMNAME_RO_HRU"|"hru"|' \
    -e 's|is_remap = true|is_remap = false|' \
    -e 's|"PARAMETER_NML"|"param.nml"|' \
    -e 's|"SEG_ID"|-9999|' \
    "${MIZUROUTE_ROOT}/route/settings/SAMPLE_toml" > "${RUN_DIR}/mizuRoute.toml"

# Force NetCDF-4/classic I/O for standalone run (avoiding Parallel NetCDF (default), which requires parallel MPI-IO and pnetcdf libraries not needed for serial standalone routing)
echo 'pio_netcdf_type = "netcdf"' >> "${RUN_DIR}/mizuRoute.toml"

echo "Running standalone mizuRoute with ${RUN_DIR}/mizuRoute.toml..."
( cd "${RUN_DIR}" && "${EXE_BIN}" "${RUN_DIR}/mizuRoute.toml" )

if ! ls "${RUN_DIR}"/standalone_smoke*.h.*.nc >/dev/null 2>&1; then
    echo "ERROR: Standalone run did not produce expected history NetCDF file in ${RUN_DIR}."
    exit 1
fi

echo "=== Standalone Build & Run Test Passed ==="
