cmake_minimum_required(VERSION 2.8)

add_subdirectory(cvc)
add_subdirectory(multi_sdf)
add_subdirectory(sign_distance_function)
add_subdirectory(xmlrpc)

# Add a special target to clean nvcc generated files.
if(CUDA_BUILD_ENABLED)
  cuda_build_clean_target()
endif(CUDA_BUILD_ENABLED)
