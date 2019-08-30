#*************************************************************************
#zdvis: Lagrangian Visualization for Vector, Tensor, and Multifield Data.
#
#Author: Zi'ang Ding
#
#Copyright (c) 2016-2018, Purdue University
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.
#**************************************************************************
cmake_minimum_required(VERSION 3.1)

set (ZD_VERSION_MAJOR 0)
set (ZD_VERSION_MINOR 3)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)

project(ZDVIS LANGUAGES C CXX)

# binary output dir

link_directories( ${LIBRARY_OUTPUT_PATH} )

if (WIN32)
  # set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/win64)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
else (WIN32)
  # set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/linux)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif(WIN32)

# dataset
set(PROJECT_DATASET_PATH "" CACHE PATH "path to the dataset")

# variables
set (OMP_MAX_THREADS "8" CACHE STRING "maximum number of threads used in OpenMP")
set (PATHLINE_MAX_COUNT "500" CACHE STRING "pre-allocate pathline vertices")

include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/nvis )

# package teem
find_package (Teem REQUIRED)
if (Teem_FOUND)
  message( "found Teem: " ${Teem_VERSION_MAJOR} "." ${Teem_VERSION_MINOR} )
  include_directories(${Teem_INCLUDE_DIRS})
else (Teem_FOUND)
  message(SEND_ERROR "Cannot find Teem")
endif (Teem_FOUND)

# package VTK
find_package (VTK REQUIRED)
if (VTK_FOUND)
  message( "found VTK: " ${VTK_MAJOR_VERSION} "." ${VTK_MINOR_VERSION} )
  include_directories(${VTK_INCLUDE_DIRS})
else (VTK_FOUND)
  message(SEND_ERROR "Cannot find VTK")
endif (VTK_FOUND)

# package Eigen
find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  message( "found Eigen3: " ${EIGEN3_VERSION_MAJOR} "." ${EIGEN3_VERSION_MINOR} )
  message(${EIGEN3_INCLUDE_DIR})
  include_directories(${EIGEN3_INCLUDE_DIR})
else (EIGEN3_FOUND)
  message(SEND_ERROR "Cannot find Eigen3")
endif (EIGEN3_FOUND)

# package NetCDF
# if (WIN32)
  find_package(netCDF)
  if (netCDF_FOUND)
    message( "found NETCDF: " ${netCDF_LIBRARIES})
    include_directories(${netCDF_INCLUDE_DIR})
  else ()
    message(SEND_ERROR "Cannot find NETCDF")
  endif ()
# else(WIN32)
#   set(netCDF_INCLUDE_DIR "" CACHE TYPE PATH)
#   set(netCDF_LIBRARIES "" CACHE TYPE FILEPATH)
# endif(WIN32)

# openmp
find_package(OpenMP)
if (OPENMP_FOUND)
  message("using OpenMP!")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-unused-result")
else (OPENMP_FOUND)
    find_package(TBB)
    if (TBB_FOUND)
        message("Using TBB")
        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_TBB")
    else(TBB_FOUND)
        message("no OpenMP support!")
    endif(TBB_FOUND)
endif(OPENMP_FOUND)


include_directories(${CMAKE_SOURCE_DIR}/utils/)
include_directories(${CMAKE_SOURCE_DIR})

macro(add_project name sourcefile)
  add_executable(${name} ${sourcefile})
  set(ExtraMacroArgs ${ARGN})
  list(LENGTH ExtraMacroArgs NumExtraMacroArgs)
  foreach(ExtraArg ${ExtraMacroArgs})
    target_sources(${name} PUBLIC ${ExtraArg})
  endforeach()

  # teem
  target_link_libraries(${name} ${Teem_LIBRARIES})
  # VTK
  target_link_libraries(${name} ${VTK_LIBRARIES})
  # netCDF
  target_link_libraries(${name} ${netCDF_LIBRARIES})
endmacro()

# configure a header file to pass some of the CMake settings to the source code
configure_file(utils/CMakeConfigure.h.in ${CMAKE_SOURCE_DIR}/utils/cmake_config.hpp)


# add_subdirectory(celltree)
add_subdirectory(Flow)
add_subdirectory(LagrangianEulerian)
add_subdirectory(Ridge)
add_subdirectory(DTI)
add_subdirectory(HARDI)
