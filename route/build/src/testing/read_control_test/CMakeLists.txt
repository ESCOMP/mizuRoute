# CMake build configuration for the read_control unit test suite
cmake_minimum_required(VERSION 3.20)

project(read_control_test LANGUAGES Fortran)

set(BUILD_TESTING OFF CACHE BOOL "Disable toml-f internal tests" FORCE)
add_subdirectory(
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../../../libraries/toml-f"
  "${CMAKE_CURRENT_BINARY_DIR}/toml-f"
)

add_executable(test.exe)

target_sources(test.exe
  PRIVATE
    "src/main.f90"
    "src/globalData.f90"
    "../../nrtype.f90"
    "../../public_var.f90"
    "../../ascii_utils.f90"
    "../../nr_utils.f90"
    "../../datetime_data.f90"
    "../../dataTypes.f90"
    "../../var_lookup.f90"
    "../../popMetadat.f90"
    "../../read_control.f90"
)

target_link_libraries(test.exe PRIVATE toml-f)

set_target_properties(test.exe PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(test.exe PRIVATE
    -g
    -Wall
    -ffree-line-length-none
    -fmax-errors=0
    -fbacktrace
    -fcheck=bounds
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(test.exe PRIVATE
    -p
    -g
    -debug
    -warn
    all
    -check
    bounds
    -noerror_limit
    -O0
    -fpe0
    -traceback
    -fltconsistency
  )
  target_link_options(test.exe PRIVATE
    -check bounds
  )
endif()

enable_testing()

add_test(
  NAME read_control_unit_test
  COMMAND $<TARGET_FILE:test.exe>
)
