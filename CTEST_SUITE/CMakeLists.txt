# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(regression_test)

include(CTest) # We include CTest which is part of CMake
include(variables.cmake)

set(ROOT ${SFEMaNS_DIR})
message (STATUS "${ROOT}")

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

#
#add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST1 ${ROOT}/REGRESSION_SUITE/TEST1/BUILD)
#add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST2 ${ROOT}/REGRESSION_SUITE/TEST2/BUILD)
#add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST3 ${ROOT}/REGRESSION_SUITE/TEST3/BUILD)
#add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST4 ${ROOT}/REGRESSION_SUITE/TEST4/BUILD)
#add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST5 ${ROOT}/REGRESSION_SUITE/TEST5/BUILD)
add_subdirectory(${ROOT}/REGRESSION_SUITE/TEST19_20 ${ROOT}/REGRESSION_SUITE/TEST19_20/BUILD)

enable_testing()
