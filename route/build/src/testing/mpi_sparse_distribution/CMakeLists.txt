cmake_minimum_required(VERSION 3.20)

project(mpi_sparse_distribution LANGUAGES Fortran)

find_package(MPI REQUIRED COMPONENTS Fortran)

#add_executable(test.exe
#  "${SRC_DIR}/main.f90"
#)

add_executable(test.exe)

target_sources(test.exe
  PRIVATE
    "src/main.f90"
    "src/globalData.f90"
    "../../mpi_utils.f90"
    "../../public_var.f90"
    "../../nrtype.f90"
)

target_link_libraries(test.exe PRIVATE MPI::MPI_Fortran)

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

#enable_testing()
#
#add_test(
#  NAME sparse_distribution_mpi_4
#  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test.exe>
#)
