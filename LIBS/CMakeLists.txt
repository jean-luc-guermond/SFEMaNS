MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Project name and version
project(sfemans_lib)

# Compile C routines and create library
enable_language(C)

file(GLOB_RECURSE sfemans_sources_C 
  ../SOURCES/WRITE_BLOCK_MM/write_block_c.c)

add_library(sfemans_lib_C SHARED ${sfemans_sources_C})

# Enable Fortran
enable_language(Fortran)

# Getting useful libraries
find_package(PkgConfig)

pkg_check_modules(PETSC PETSc)
pkg_check_modules(FFTW fftw3)

# Build library sfemans_lib_F
file(GLOB_RECURSE sfemans_sources_F
	../SOURCES/FEMSUB/*.F90
	../SOURCES/FEMSUB_OBJET/*.F90
	../SOURCES/GAUSS_POINTS_OBJET/*.F90
	../SOURCES/FEMSUB_PARALLEL/*.F90
	../SOURCES/FEMSUB_SFEMaNS/*.F90
	../SOURCES/MHD_DATA_TEST_CONV_PETSC/*.F90
	../SOURCES/WRITE_BLOCK_MM/*.F90
  )
add_library(sfemans_lib_F SHARED ${sfemans_sources_F})
target_include_directories(sfemans_lib_F PUBLIC
    ${FFTW_INCLUDEDIR} ${FFTW_INCLUDE_DIRS}
    ${PETSC_INCLUDEDIR} ${PETSC_INCLUDE_DIRS}
  )
target_link_libraries(sfemans_lib_F sfemans_lib_C) 
