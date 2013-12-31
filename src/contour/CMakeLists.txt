cmake_minimum_required (VERSION 2.6)

file(GLOB SOURCE_FILES *.cpp *.c)
file(GLOB INCLUDE_FILES ../../inc/Contour/*.[Hh])

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

include_directories(
  ../../inc
)

add_definitions (-fpermissive -D_BOOL -DSP2 -DUSEDICT -DNDEBUG)
if(WIN32)
  add_definitions(-DWIN32)
endif(WIN32)
if(APPLE)
  add_definitions(-DMACOS_X)
endif(APPLE)

add_library(contour STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

#target_link_libraries(contour Utility)
