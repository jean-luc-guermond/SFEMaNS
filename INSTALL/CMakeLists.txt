MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Includes various variables depending on native
include(variables.cmake)

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Project name and version
project(debug_sfemans)

# Compile C routines and create library
enable_language(C)

file(GLOB_RECURSE sfemans_sources_C 
  ${SFEMaNS_DIR}/WRITE_BLOCK_MM/write_block_c.c
  )

add_library(sfemans_lib_C SHARED ${sfemans_sources_C})


# Enable Fortran
enable_language(Fortran)
add_compile_options(${COMPILE_OPTIONS})


# Getting useful libraries
find_package(MPI)
find_package(ZLIB)
find_package(PkgConfig)

pkg_check_modules(PETSC PETSc)
pkg_check_modules(FFTW fftw3)

# Searching all petsc linked librariesi
execute_process (COMMAND pkg-config PETSc --libs --static OUTPUT_VARIABLE PETSC_PRIVATE_LIBRARIES)
string(STRIP "${PETSC_PRIVATE_LIBRARIES}" PETSC_PRIVATE_LIBRARIES)

# Setting for different builds
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE " -O3 " CACHE STRING
  " Flags used by the compiler during release builds. ")
set(CMAKE_Fortran_FLAGS_NATIVE " -O4 -march=native -mtune=native -Ofast " CACHE STRING
  " Flags used by the compiler during release builds. ")

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Native"
  )



# Since the SFEMaNS directory *should* not change structure, i.e.
# files are not added/removed, we can use the file()
# command to populate a variable with the list of needed files

file(GLOB_RECURSE sfemans_sources 
  ${SFEMaNS_DIR}/FEMSUB/*.f90
  ${SFEMaNS_DIR}/FEMSUB_OBJET/*.f90
  ${SFEMaNS_DIR}/GAUSS_POINTS_OBJET/*.f90
  ${SFEMaNS_DIR}/FEMSUB_PARALLEL/*.f90
  ${SFEMaNS_DIR}/FEMSUB_SFEMaNS/*.f90
  ${SFEMaNS_DIR}/MHD_DATA_TEST_CONV_PETSC/*.f90
  ${SFEMaNS_DIR}/WRITE_BLOCK_MM/*.f90
  )

file(GLOB_RECURSE local_sources ${SFEMaNS_DIR}/TEMPLATE/*.f90)

# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../EXECUTABLE)

# Add executable
add_executable(a.exe ${sfemans_sources} ${local_sources})

# Include and link library
target_include_directories(a.exe PUBLIC ${FFTW_INCLUDEDIR} ${PETSC_INCLUDEDIR})
target_link_libraries(a.exe PUBLIC sfemans_lib_C ${PETSC_LINK_LIBRARIES} ${PETSC_PRIVATE_LIBRARIES} ${FFTW_LINK_LIBRARIES} ${ZLIB_LIBRARIES} ${MPI_Fortran_LIBRARIES})
