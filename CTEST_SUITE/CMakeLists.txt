# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

include(variables.cmake)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(regression_test)

include(CTest) # We include CTest which is part of CMake

set(ROOT ${SFEMaNS_DIR})
message (STATUS "${ROOT}")


set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug_bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE ${release_bounds} CACHE STRING
	  " Flags used by the compiler during release builds. ")
set(CMAKE_Fortran_FLAGS_NATIVE ${native_bounds} CACHE STRING
	    " Flags used by the compiler during native builds. ")

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug")

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release")

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Native")

# Getting useful libraries
find_package(MPI)
find_package(ZLIB)
find_package(PkgConfig)

pkg_check_modules(PETSC PETSc)
pkg_check_modules(FFTW fftw3)

# Adding directory to build lib_f
add_subdirectory(${SFEMaNS_DIR}/LIBS ${SFEMaNS_DIR}/LIBS)

add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST1 ${ROOT}/REGRESSION_SUITE/TEST1/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST2 ${ROOT}/REGRESSION_SUITE/TEST2/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST3 ${ROOT}/REGRESSION_SUITE/TEST3/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST4 ${ROOT}/REGRESSION_SUITE/TEST4/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST5 ${ROOT}/REGRESSION_SUITE/TEST5/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST6 ${ROOT}/REGRESSION_SUITE/TEST6/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST7 ${ROOT}/REGRESSION_SUITE/TEST7/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST8 ${ROOT}/REGRESSION_SUITE/TEST8/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST9 ${ROOT}/REGRESSION_SUITE/TEST9/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST10 ${ROOT}/REGRESSION_SUITE/TEST10/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST11 ${ROOT}/REGRESSION_SUITE/TEST11/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST12 ${ROOT}/REGRESSION_SUITE/TEST12/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST13 ${ROOT}/REGRESSION_SUITE/TEST13/BUILD)
#TEST 14 Disabled (=always passes) to avoid installing arpack
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST14 ${ROOT}/REGRESSION_SUITE/TEST14/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST15 ${ROOT}/REGRESSION_SUITE/TEST15/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST16 ${ROOT}/REGRESSION_SUITE/TEST16/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST17 ${ROOT}/REGRESSION_SUITE/TEST17/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST18 ${ROOT}/REGRESSION_SUITE/TEST18/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST19_20 ${ROOT}/REGRESSION_SUITE/TEST19_20/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST21 ${ROOT}/REGRESSION_SUITE/TEST21/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST22 ${ROOT}/REGRESSION_SUITE/TEST22/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST23 ${ROOT}/REGRESSION_SUITE/TEST23/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST24 ${ROOT}/REGRESSION_SUITE/TEST24/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST25 ${ROOT}/REGRESSION_SUITE/TEST25/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST26 ${ROOT}/REGRESSION_SUITE/TEST26/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST27 ${ROOT}/REGRESSION_SUITE/TEST27/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST28 ${ROOT}/REGRESSION_SUITE/TEST28/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST29 ${ROOT}/REGRESSION_SUITE/TEST29/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST30 ${ROOT}/REGRESSION_SUITE/TEST30/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST31_32 ${ROOT}/REGRESSION_SUITE/TEST31_32/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST33 ${ROOT}/REGRESSION_SUITE/TEST33/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST34 ${ROOT}/REGRESSION_SUITE/TEST34/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST35 ${ROOT}/REGRESSION_SUITE/TEST35/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST36 ${ROOT}/REGRESSION_SUITE/TEST36/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST37 ${ROOT}/REGRESSION_SUITE/TEST37/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST38 ${ROOT}/REGRESSION_SUITE/TEST38/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST39 ${ROOT}/REGRESSION_SUITE/TEST39/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST40 ${ROOT}/REGRESSION_SUITE/TEST40/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST41 ${ROOT}/REGRESSION_SUITE/TEST41/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST42 ${ROOT}/REGRESSION_SUITE/TEST42/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST43_44 ${ROOT}/REGRESSION_SUITE/TEST43_44/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST45 ${ROOT}/REGRESSION_SUITE/TEST45/BUILD)

enable_testing()
