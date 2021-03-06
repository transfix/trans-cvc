cmake_minimum_required(VERSION 2.8)
project(trans-cvc)

# CMake Modules
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/CMake/cuda
  ${CMAKE_MODULE_PATH})

# Build output 
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# The version number.
set(transcvc_VERSION_MAJOR 1)
set(transcvc_VERSION_MINOR 0)

# The full version string
set(transcvc_VERSION ${transcvc_VERSION_MAJOR}.${transcvc_VERSION_MINOR})
message(STATUS "trans-cvc version: ${transcvc_VERSION}")

# include macros for setting up some library dependencies
#include( SetupQt )
include(SetupFFTW)
include(SetupBoost)
include(SetupGSL)

option(DISABLE_CGAL "This disables CGAL even if it is available on the build system." OFF)
mark_as_advanced(DISABLE_CGAL)
if (NOT DISABLE_CGAL)
  include( SetupCGAL )
endif (NOT DISABLE_CGAL)

# 05/11/2012 - transfix
# Possibly prevent crashing on OSX Lion
# http://groups.google.com/group/boost-list/browse_thread/thread/1c5a614a9e8c9808
if(APPLE)
  add_definitions(-D_GLIBCXX_FULLY_DYNAMIC_STRING)
endif(APPLE)

# experimental
option(USING_STANDARD_INSTALL_LOCATION "If on, it will not put everything under $prefix/trans-cvc" ON)
mark_as_advanced(USING_STANDARD_INSTALL_LOCATION)
if(USING_STANDARD_INSTALL_LOCATION)
  set(ADDED_PREFIX ".")
else(USING_STANDARD_INSTALL_LOCATION)
  set(ADDED_PREFIX "trans-cvc")
endif(USING_STANDARD_INSTALL_LOCATION)

add_subdirectory (src)
