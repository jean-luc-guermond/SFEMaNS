MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Includes various variables depending on native
include(variables.cmake)

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Project name and version
project(debug_sfemans)

# Enable Fortran
enable_language(Fortran)

# Getting useful libraries
find_package(MPI)
find_package(ZLIB)
find_package(PkgConfig)

pkg_check_modules(PETSC PETSc)
pkg_check_modules(FFTW fftw3)

# Searching all petsc linked libraries
execute_process (COMMAND pkg-config PETSc --libs --static OUTPUT_VARIABLE PETSC_PRIVATE_LIBRARIES)
string(STRIP "${PETSC_PRIVATE_LIBRARIES}" PETSC_PRIVATE_LIBRARIES)

# Create local sources from DEBUG
file(GLOB_RECURSE template_sources ${SFEMaNS_DIR}/TEMPLATE/*.F90)

# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./EXECUTABLE)

# Add executable
add_executable(a.exe ${template_sources})

# Include and link library
target_link_libraries(a.exe PUBLIC sfemans_lib_F
  ${PETSC_LINK_LIBRARIES}
  ${PETSC_PRIVATE_LIBRARIES}
  ${FFTW_LINK_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${MPI_Fortran_LIBRARIES}
  ${ADDITIONAL_LINKS})

target_include_directories(a.exe PUBLIC ${SFEMaNS_DIR}/LIBS)
# REGRESSION TESTS

# Create local sources from TEMPLATE
file(GLOB_RECURSE regression_sources
  ./condlim.F90
  ./read_user_data.F90
  ${SFEMaNS_DIR}/SOURCES/MESH_INTERPOLATION/*.F90)

# Add executable
add_executable(a_mesh_interpolation.exe ${regression_sources})

# Include and link library
target_link_libraries(a_mesh_interpolation.exe PUBLIC sfemans_lib_F
  ${PETSC_LINK_LIBRARIES}
  ${PETSC_PRIVATE_LIBRARIES}
  ${FFTW_LINK_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${MPI_Fortran_LIBRARIES}
  ${ADDITIONAL_LINKS})

target_include_directories(a_mesh_interpolation.exe PUBLIC ${SFEMaNS_DIR}/LIBS)

ADD_CUSTOM_TARGET(regression ALL
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target a.exe -j20
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target a_mesh_interpolation.exe -j20
  COMMAND ./debug_SFEMaNS_template ${SFEMaNS_DIR} ./EXECUTABLE/a.exe ./EXECUTABLE/a_mesh_interpolation.exe ${RUN_PRE_PROC} ${PROC_CALL} ${RUN_POST_PROC}
)

# Adding directory to build lib_f
add_subdirectory(${SFEMaNS_DIR}/LIBS ${SFEMaNS_DIR}/LIBS)
