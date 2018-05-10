cmake_minimum_required(VERSION 2.8)

add_subdirectory(cvc)

option(CVC_USING_XMLRPC "If this is enabled, cvc::state will support state sharing with other processes on the network" OFF)
mark_as_advanced(CVC_USING_XMLRPC)
if(CVC_USING_XMLRPC)
  add_subdirectory(xmlrpc)
endif()

# Add a special target to clean nvcc generated files.
if(CUDA_BUILD_ENABLED)
  cuda_build_clean_target()
endif(CUDA_BUILD_ENABLED)
