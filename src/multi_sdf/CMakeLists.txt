cmake_minimum_required(VERSION 2.6)

set(SOURCE_FILES
  dt.cpp
  init.cpp
  tcocone.cpp
  robust_cc.cpp
  kdtree.cpp
  mesh_io.cpp
  op.cpp
  sdf.cpp
  util.cpp
  rcocone.cpp
  priorityqueue.cpp
  matrix.cpp
)
file(GLOB INCLUDE_FILES ../../inc/multi_sdf/*.[Hh])

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

include_directories(
  ../../inc
)

add_library(multi_sdf STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

find_package(CGAL)
if(CGAL_FOUND)
  include(${CGAL_USE_FILE})
  # need the following flags in case CGAL has some special compiler needs for this compiler
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")
  set(LIBS ${CGAL_LIBRARIES})
else(CGAL_FOUND)
  message(SEND_ERROR "multi_sdf requires CGAL!")
endif(CGAL_FOUND)

target_link_libraries(multi_sdf
  ${LIBS}
)

option(BUILD_MULTISDF_STANDALONE "Build standalone multi_sdf program." OFF)
mark_as_advanced(BUILD_MULTISDF_STANDALONE)
if(BUILD_MULTISDF_STANDALONE)
  source_group("Source Files" FILES ${SOURCE_FILES} main.cpp)
  add_executable(multiSDF main.cpp)
  target_link_libraries(multiSDF multi_sdf)    
endif(BUILD_MULTISDF_STANDALONE)
