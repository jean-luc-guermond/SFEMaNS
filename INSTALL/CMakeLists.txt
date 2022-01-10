MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Project name and version
project(debug_sfemans)

# Compile C routines and create library
enable_language(C)
file(GLOB_RECURSE sfemans_sources_C 
  $ENV{SFEMaNS_DIR}/WRITE_BLOCK_MM/write_block_c.c
  )

add_library(sfemans_lib_C SHARED ${sfemans_sources_C})


# Enable Fortran
enable_language(Fortran)
add_compile_options(-x f95-cpp-input)

find_package(MPI)
find_package(HDF5)
find_package(ZLIB)
find_package(PkgConfig)

pkg_check_modules(PETSC PETSc)
pkg_check_modules(FFTW fftw3)

include_directories(${PETSC_INCLUDE_DIRS})
#set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")     

execute_process (COMMAND pkg-config PETSc --libs --static OUTPUT_VARIABLE PETSC_PRIVATE_LIBRARIES)
string(STRIP "-L /usr/lib/x86_64-linux-gnu/hdf5/serial ${PETSC_PRIVATE_LIBRARIES}" PETSC_PRIVATE_LIBRARIES)

#include_directories(${FFTW_INCLUDES})
message(${FFTW_INCLUDES} " FFTW INCLUDES")
include_directories(/usr/include)
message(${FFTW_LINK_LIBRARIES} " FFTW LIBRARIES")



#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
#set(CMAKE_Fortran_FLAGS_RELEASE " -O3 " CACHE STRING
#  " Flags used by the compiler during release builds. ")
#set(CMAKE_Fortran_FLAGS_NATIVE " -O4 -march=native -mtune=native -Ofast " CACHE STRING
#  " Flags used by the compiler during release builds. ")

# Define runtime output directory
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../EXECUTABLE)

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
  $ENV{SFEMaNS_DIR}/FEMSUB/*.f90
  $ENV{SFEMaNS_DIR}/FEMSUB_OBJET/*.f90
  $ENV{SFEMaNS_DIR}/GAUSS_POINTS_OBJET/*.f90
  $ENV{SFEMaNS_DIR}/FEMSUB_PARALLEL/*.f90
  $ENV{SFEMaNS_DIR}/FEMSUB_SFEMaNS/*.f90
  $ENV{SFEMaNS_DIR}/MHD_DATA_TEST_CONV_PETSC/*.f90
  $ENV{SFEMaNS_DIR}/WRITE_BLOCK_MM/*.f90
  )

file(GLOB_RECURSE local_sources ${CMAKE_SOURCE_DIR}/LOCAL_SOURCES/*.f90)


# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../EXECUTABLE)

# add executable
add_executable(a.exe
  ${sfemans_sources}
  ${local_sources}
  )

# include and link library
target_link_libraries(a.exe  PUBLIC sfemans_lib_C ${PETSC_LINK_LIBRARIES} ${PETSC_PRIVATE_LIBRARIES} ${FFTW_LINK_LIBRARIES}
  ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${HDF5_LIBRARIES} -lmetis ${ZLIB_LIBRARIES} -lz)
