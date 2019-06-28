###################################################################################
#
# CubeZ
#
# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
#
#-*- CMAKE -*-#
# Process this file with cmake to produce a configure script.

# Following parameters must be given at cmake command line
# INSTALL_DIR, CXX, CC, F90

# Build options, see Readme.md for details.
#
# -D INSTALL_DIR=install_directory
#
# -D with_MPI=(yes|no)
#
# -D real_type=(float|double)
#
# -D enable_OPENMP=(yes|no)
#
# -D with_PM=installed_directory or OFF
#
# -D with_CBR=installed_directory or OFF
#
# -D with_PAPI=installed_directory or OFF
#
# -D with_SIMD=(OFF|256|512)
#
# -D with_Ftrace=(off|on)
#
# -D with_ACC={off|Pascal|Volta}
#

cmake_minimum_required(VERSION 2.6)


#######
# Build options
#######

option(with_MPI "Enable MPI" "ON")
option(enable_OPENMP "Enable OpenMP" "ON")
option(real_type "Type of floating point" "OFF")
option(with_PAPI "Enable PAPI" "OFF")
option(with_SIMD "SIMD Width" "OFF")
option(trace "with_Ftrace Option" "OFF")
option(with_ACC "OpenACC Option" "OFF")

#######
# Project setting
#######

project(CubeZ C CXX Fortran)

if(INSTALL_DIR STREQUAL "")
  set(INSTALL_DIR "/usr/local/CubeZ")
endif()

# project()でCMAKE_INSTALL_PREFIX=/usr/localとなるので、再定義
set (CMAKE_INSTALL_PREFIX "${INSTALL_DIR}")

set(CMAKE_BUILD_TYPE Release)

#######
# Policy setting
#######

# $ cmake --help-policy CMP0012
cmake_policy(SET CMP0012 NEW)


#######
set(PROJECT_VERSION "1.0.3")
set(LIB_REVISION "20190627_2300")
#######



#######
# Cmake module
#######

set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileOptionSelector.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/cz_package.cmake)



#######
# Options
#######

AddOptimizeOption()

# Fortran free form
FreeForm()


# OpennMP
checkOpenMP()


# Real type
precision()

# Trace iff Aurora
if(TARGET_ARCH STREQUAL "NEC_Aurora_VE")
  if(with_Ftrace)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrace")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ftrace")
  endif()
endif()


# Serial flag
if(NOT with_MPI)
  SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -DDISABLE_MPI")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_MPI")
endif()



if(USE_F_TCS STREQUAL "YES")
else()
  if(with_MPI)
    include(FindMPI)

    if(MPI_CXX_FOUND)
      message("MPI-CXX-CMP  : ${MPI_CXX_COMPILER}")
      message("MPI-CXX-FLAGS: ${MPI_CXX_COMPILE_FLAGS}")
      message("MPI-CXX-INC  : ${MPI_CXX_INCLUDE_PATH}")
      message("MPI-LNK-FLAGS: ${MPI_CXX_LINK_FLAGS}")
      message("MPI-LIBS     : ${MPI_CXX_LIBRARIES}")
    endif(MPI_CXX_FOUND)
  endif()
endif()

message(" ")
message("PROJECT : CubeZ " ${PROJECT_VERSION} " : " ${LIB_REVISION})
message(" ")



#######
# find Libraries to depend on
#######


# PMlib
if(with_PM STREQUAL "OFF")
  SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -DDISABLE_PMLIB")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_PMLIB")
elseIF(with_PM)
  SET(PM_DIR "${with_PM}")
  SET(PM_INC "${PM_DIR}/include")
  SET(PM_LIB "${PM_DIR}/lib")
ELSE()
  MESSAGE("Error: can not find PMlib.")
ENDIF()


# PAPI
if(with_PAPI STREQUAL "OFF")
elseif(with_PAPI)
  SET(PAPI_DIR "${with_PAPI}")
  SET(PAPI_INC "${PAPI_DIR}/include")
  SET(PAPI_LIB "${PAPI_DIR}/lib")
else()
  MESSAGE("Error: can not find PAPI")
endif()


# CBrick
IF(with_CBR STREQUAL "OFF")
elseIF(with_CBR)
  SET(CBR_DIR "${with_CBR}")
  SET(CBR_INC "${CBR_DIR}/include")
  SET(CBR_LIB "${CBR_DIR}/lib")
ELSE()
  MESSAGE("Error: can not find CBrick.")
ENDIF()


####
# OpenACC
####

if(with_ACC)
  if (with_ACC STREQUAL "Pascal")
    SET(ACC "-ta=tesla,cc60")
  elseif (with_ACC STREQUAL "Volta")
    SET(ACC "-ta=tesla,cc70")
  endif()
endif()




#######
# Display options
#######

message(" ")
message("${PROJECT_NAME} Version : ${PROJECT_VERSION}")
message(" ")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "CMAKE_CXX_FLAGS       : "      ${CMAKE_CXX_FLAGS})
message( STATUS "CMAKE_Fortran_COMPILER: "      ${CMAKE_Fortran_COMPILER})
message( STATUS "CMAKE_Fortran_FLAGS   : "      ${CMAKE_Fortran_FLAGS})

message( STATUS "MPI support           : "      ${with_MPI})
message( STATUS "OpenMP support        : "      ${enable_OPENMP})
message( STATUS "Type of floating point: "      ${real_type})
message( STATUS "PMlib support         : "      ${with_PM})
message( STATUS "PAPI support          : "      ${with_PAPI})
message( STATUS "CBrick support        : "      ${with_CBR})
message( STATUS "SIMD WIDTH            : "      ${with_SIMD})
message( STATUS "OpenACC support       : "      ${with_ACC} " " ${ACC})
message( STATUS "Aurora Ftrace         : "      ${with_Ftrace})
message(" ")


#######
# Sub Directories
#######

add_subdirectory(src)
#add_subdirectory(example)



#######
# configure files
#######

configure_file( src/cz_cpp/czVersion.h.in ${PROJECT_BINARY_DIR}/src/cz_cpp/czVersion.h @ONLY)



#######
# Install files
#######

install(FILES  ChangeLog.md
               Readme.md
               License.txt
               ${PROJECT_BINARY_DIR}/CMakeCache.txt
        DESTINATION share)


#######
# Example
#######

#add_subdirectory(example/tri)


#######
# Installer
#######

cz_package()
