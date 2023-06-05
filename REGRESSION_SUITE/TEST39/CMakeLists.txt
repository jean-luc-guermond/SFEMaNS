MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Includes various variables depending on native
#include(variables.cmake)

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Project name and version
project(sfemans)

# Enable Fortran
enable_language(Fortran)

# Setting for different builds
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug_bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE ${release_bounds} CACHE STRING
  " Flags used by the compiler during release builds. ")
set(CMAKE_Fortran_FLAGS_NATIVE ${native_bounds} CACHE STRING
  " Flags used by the compiler during native builds. ")

IF (TARGET debug)
else()
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all -j20
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all -j20
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all -j20
  COMMENT "Switch CMAKE_BUILD_TYPE to Native"
  )
endif()

# Create local sources from local dir
file(GLOB_RECURSE local_sources ${SFEMaNS_DIR}/TEMPLATE/*.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/*.F90)
list(REMOVE_ITEM local_sources ${SFEMaNS_DIR}/TEMPLATE/condlim.F90)

# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../EXECUTABLE)

# Add executable
set(exe "test39.exe")
add_executable(${exe} ${local_sources})

# Include and link library
target_link_libraries(${exe} PUBLIC sfemans_lib_F
  ${PETSC_LINK_LIBRARIES}
  ${PETSC_PRIVATE_LIBRARIES}
  ${FFTW_LINK_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${MPI_Fortran_LIBRARIES}
  ${ADDITIONAL_LINKS}
  )

target_include_directories(${exe} PUBLIC ${SFEMaNS_DIR}/LIBS)

set(passRegex "123")

set (nm "test39")
add_test(NAME ${nm}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../REGRESSION_TESTS
  COMMAND /bin/bash/ job.sh ${RUN_PRE_PROC} ${PROC_CALL} ${exe} ${RUN_POST_PROC})
set_property (TEST ${nm} PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
