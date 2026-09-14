cmake_minimum_required(VERSION 3.20)

project(netcdf_test LANGUAGES Fortran)

find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDFF REQUIRED IMPORTED_TARGET netcdf-fortran)


add_executable(test.exe)

target_sources(test.exe
  PRIVATE
    "src/example_netcdf.f90"
    "src/test_netcdf.f90"
)

target_link_libraries(test.exe PRIVATE PkgConfig::NETCDFF)

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
