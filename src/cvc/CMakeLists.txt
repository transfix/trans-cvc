cmake_minimum_required(VERSION 2.6)

file(GLOB INCLUDE_FILES ../../inc/cvc/*.[Hh])
set(SOURCE_FILES
  app.cpp
  state.cpp
)

include_directories(
  ../../inc
  ${CMAKE_BINARY_DIR}/inc
)

set(CVC_NAMESPACE cvc)
mark_as_advanced(CVC_NAMESPACE)
add_definitions(-DCVC_NAMESPACE=${CVC_NAMESPACE})

# The libcvc version number.
set(CVC_VERSION_MAJOR 1)
set(CVC_VERSION_MINOR 0)

# The full version string
set(CVC_VERSION ${CVC_VERSION_MAJOR}.${CVC_VERSION_MINOR})
message(STATUS "libcvc version: ${CVC_VERSION}")

add_definitions(-DCVC_VERSION_STRING=\"${CVC_VERSION}\")

set(CVC_LINK_LIBS "")
option(CVC_USING_XMLRPC "If this is enabled, cvc::state will support state sharing with other processes on the network" OFF)
if(CVC_USING_XMLRPC)
  set(CVC_LINK_LIBS xmlrpc)
  add_definitions(-DUSING_XMLRPC)
endif(CVC_USING_XMLRPC)

option(CVC_USING_HDF5 "HDF5 Support. Needed for *.cvc format but requires linking to libhdf5." ON)
mark_as_advanced(CVC_USING_HDF5)
if(CVC_USING_HDF5)
  find_package(HDF5 COMPONENTS C CXX)
  if(HDF5_FOUND)
    message("HDF5 found! Enabling *.cvc file support")

    add_definitions(-DCVC_USING_HDF5 ${HDF5_DEFINITIONS})
    set(INCLUDE_FILES
       ${INCLUDE_FILES}
       ../../inc/cvc/hdf5_utils.h
    )
    set(SOURCE_FILES
       ${SOURCE_FILES}
       hdf5_utils.cpp
    )

    include_directories(${HDF5_INCLUDE_DIR})

    set(CVC_LINK_LIBS ${CVC_LINK_LIBS} ${HDF5_LIBRARIES})
  else(HDF5_FOUND)
    message("HDF5 not found! *.cvc file support disabled")
    set(CVC_HDF5_DISABLED TRUE)
  endif(HDF5_FOUND)
else(CVC_USING_HDF5)
  set(CVC_HDF5_DISABLED TRUE)
endif(CVC_USING_HDF5)

option(CVC_LOG4CPLUS_DEFAULT "Use log4cplus output by default" OFF)
mark_as_advanced(CVC_LOG4CPLUS_DEFAULT)
if(CVC_LOG4CPLUS_DEFAULT)
  find_package(Log4cplus)
  if(LOG4CPLUS_FOUND)
    add_definitions(-DUSING_LOG4CPLUS_DEFAULT)
    include_directories(${LOG4CPLUS_INCLUDE_DIR})
    set(CVC_LINK_LIBS ${CVC_LINK_LIBS} ${LOG4CPLUS_LIBRARIES})
  else(LOG4CPLUS_FOUND)
    message("Log4cplus not found! Log4cplus based logging disabled")
  endif(LOG4CPLUS_FOUND)
endif(CVC_LOG4CPLUS_DEFAULT)

option(LOG4CPLUS_DISABLE_TRACE "Disable logging trace" OFF)
mark_as_advanced(LOG4CPLUS_DISABLE_TRACE)

option(CVC_APP_XML_PROPERTY_MAP "Read/Write XML for property map info.  If not set, uses boost INFO property tree format." OFF)
if(CVC_APP_XML_PROPERTY_MAP)
  add_definitions(-DCVC_APP_XML_PROPERTY_MAP)
endif(CVC_APP_XML_PROPERTY_MAP)
mark_as_advanced(CVC_APP_XML_PROPERTY_MAP)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/../../inc/cvc/config.h.cmake ${CMAKE_BINARY_DIR}/inc/cvc/config.h )

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

add_library(cvc SHARED ${SOURCE_FILES} ${INCLUDE_FILES})

SetupBoost(cvc)

target_link_libraries(cvc ${CVC_LINK_LIBS})
