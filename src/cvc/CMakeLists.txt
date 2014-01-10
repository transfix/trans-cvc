cmake_minimum_required(VERSION 2.6)

if(WIN32)
  add_definitions(-D__WINDOWS__)
endif(WIN32)

file(GLOB INCLUDE_FILES ../../inc/cvc/*.[Hh])
set(SOURCE_FILES
  algorithm.cpp
  anisotropic_diffusion.cpp
  app.cpp
  bilateral_filter.cpp
  contrast_enhancement.cpp
  gdtv_filter.cpp
  geometry.cpp
  mrc_io.cpp
  null_io.cpp
  rawiv_io.cpp
  rawv_io.cpp
  smoothing.cpp
  spider_io.cpp
  state.cpp
  utility.cpp
  voxels.cpp
  volume.cpp
  volume_file_io.cpp
  volume_file_info.cpp
  vtk_io.cpp
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

set(CVC_LINK_LIBS ""
#    sign_distance_function 
#    fast_contouring
#    contour
)

option(CVC_USING_XMLRPC "If this is enabled, cvc::state will support state sharing with other processes on the network" ON)
mark_as_advanced(CVC_USING_XMLRPC)
if(CVC_USING_XMLRPC)
  set(CVC_LINK_LIBS ${CVC_LINK_LIBS} xmlrpc)
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
       hdf5_io.cpp
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

option(CVC_USING_IMOD_MRC "Add support for using imod's MRC loading routines." OFF)
mark_as_advanced(CVC_USING_IMOD_MRC)
if(CVC_USING_IMOD_MRC)
  add_definitions(-DNOTIFFLIBS -DCVC_USING_IMOD_MRC)
  include_directories(libiimod)
  file(GLOB IMOD_MRC_LIBIIMOD_HEADERS libiimod/*.h)
  set(INCLUDE_FILES
     ${INCLUDE_FILES}
     ${IMOD_MRC_LIBIIMOD_HEADERS}
  )
  set(SOURCE_FILES
     ${SOURCE_FILES}
     libiimod/b3dutil.c
     libiimod/diffusion.c
     libiimod/iilikemrc.c
     libiimod/iimage.c
     libiimod/iimrc.c
     libiimod/iitif.c
     libiimod/ilist.c
     libiimod/mrcfiles.c
     libiimod/mrcsec.c
     libiimod/mrcslice.c
     libiimod/plist.c
     libiimod/sliceproc.c
     libiimod/tiffstub.c
     libiimod/islice.c
     libiimod/parallelwrite.c
  )
endif(CVC_USING_IMOD_MRC)

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

option(CVC_GEOMETRY_ENABLE_BUNNY "Enable built-in bunny geometry test tri mesh" ON)
mark_as_advanced(CVC_GEOMETRY_ENABLE_BUNNY)
if(CVC_GEOMETRY_ENABLE_BUNNY)
  add_definitions(-DCVC_GEOMETRY_ENABLE_BUNNY)
endif(CVC_GEOMETRY_ENABLE_BUNNY)

option(CVC_GEOMETRY_CORRECT_INDEX_START "adjust indices if they start from 1 rather than 0." OFF)
mark_as_advanced(CVC_GEOMETRY_CORRECT_INDEX_START)
if(CVC_GEOMETRY_ENABLE_BUNNY)
  add_definitions(-DCVC_GEOMETRY_ENABLE_BUNNY)
endif(CVC_GEOMETRY_ENABLE_BUNNY)

if(NOT DISABLE_CGAL)
  find_package(CGAL)
  if(CGAL_FOUND)
    include(${CGAL_USE_FILE})
    # need the following flags in case CGAL has some special compiler needs for this compiler
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")
    set(CVC_LINK_LIBS ${CVC_LINK_LIBS} ${CGAL_LIBRARIES})
    add_definitions(-DCVC_GEOMETRY_ENABLE_PROJECT)
    set(INCLUDE_FILES 
        ${INCLUDE_FILES} 
	../../inc/cvc/project_verts.h)
    set(CVC_LINK_LIBS ${CVC_LINK_LIBS} ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})
  endif(CGAL_FOUND)
endif(NOT DISABLE_CGAL)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/../../inc/cvc/config.h.cmake ${CMAKE_BINARY_DIR}/inc/cvc/config.h )

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${INCLUDE_FILES})

option(CVC_ENABLE_MESHER "enable isosurfacing/meshing capability" ON)
mark_as_advanced(CVC_ENABLE_MESHER)
if(CVC_ENABLE_MESHER)
  add_definitions(-DCVC_ENABLE_MESHER -D_BOOL -DSP2 -DUSEDICT -DNDEBUG)
  include_directories(cvc-mesher/contour
                      cvc-mesher/FastContouring
                      cvc-mesher/VolMagick
                      cvc-mesher/LBIE)
  set(SOURCE_FILES
      ${SOURCE_FILES}
      cvc-mesher/contour/squeue.cpp
      cvc-mesher/contour/rangesweep.cpp
      cvc-mesher/contour/dirseeds.cpp
      cvc-mesher/contour/contour.cpp
      cvc-mesher/contour/seeddirreg3.cpp
      cvc-mesher/contour/contour3d.cpp
      cvc-mesher/contour/conplotreg3.cpp
      cvc-mesher/contour/dict.c
      cvc-mesher/contour/edgehash.cpp
      cvc-mesher/contour/regprop2.cpp
      cvc-mesher/contour/spectrumtest.cpp
      cvc-mesher/contour/ipqueue.cpp
      cvc-mesher/contour/seedcells.cpp
      cvc-mesher/contour/conplot3d.cpp
      cvc-mesher/contour/dataslc.cpp
      cvc-mesher/contour/cubes.cpp
      cvc-mesher/contour/conplot2d.cpp
      cvc-mesher/contour/shelf.cpp
      cvc-mesher/contour/iqueue.cpp
      cvc-mesher/contour/datavol.cpp
      cvc-mesher/contour/cellqueue.cpp
      cvc-mesher/contour/bin.cpp
      cvc-mesher/contour/datareg2.cpp
      cvc-mesher/contour/bucketsearch.cpp
      cvc-mesher/contour/seedall.cpp
      cvc-mesher/contour/range.cpp
      cvc-mesher/contour/hash.cpp
      cvc-mesher/contour/conplotreg2.cpp
      cvc-mesher/contour/data.cpp
      cvc-mesher/contour/inttree.cpp
      cvc-mesher/contour/seedchkr3.cpp
      cvc-mesher/contour/conplot.cpp
      cvc-mesher/contour/queue.cpp
      #cvc-mesher/contour/main.cpp
      cvc-mesher/contour/regprop.cpp
      cvc-mesher/contour/respprop2.cpp
      cvc-mesher/contour/rangeprop.cpp
      cvc-mesher/contour/segtree.cpp
      cvc-mesher/contour/conplot_p.cpp
      cvc-mesher/contour/contour2d.cpp
      cvc-mesher/contour/datareg3.cpp
      cvc-mesher/contour/dirseedsreg2.cpp
      cvc-mesher/contour/seedchkr2.cpp
      cvc-mesher/FastContouring/Ray.cpp
      cvc-mesher/FastContouring/FastContouring.cpp
      cvc-mesher/FastContouring/Quaternion.cpp
      cvc-mesher/FastContouring/Tuple.cpp
      cvc-mesher/FastContouring/Vector.cpp
      cvc-mesher/FastContouring/ContourGeometry.cpp
      cvc-mesher/FastContouring/Matrix.cpp
      cvc-mesher/FastContouring/MarchingCubesBuffers.cpp
      cvc-mesher/VolMagick/RawV.cpp
      cvc-mesher/VolMagick/MRC.cpp
      cvc-mesher/VolMagick/VolumeCache.cpp
      cvc-mesher/VolMagick/ContrastEnhancement.cpp
      cvc-mesher/VolMagick/AnisotropicDiffusion.cpp
      cvc-mesher/VolMagick/BilateralFilter.cpp
      cvc-mesher/VolMagick/INR.cpp
      cvc-mesher/VolMagick/RawIV.cpp
      cvc-mesher/VolMagick/VolMagick.cpp
      cvc-mesher/VolMagick/Spider.cpp
      cvc-mesher/VolMagick/main.cpp
      cvc-mesher/LBIE/normalspline.cpp
      cvc-mesher/LBIE/e_face.cpp
      cvc-mesher/LBIE/pcio.cpp
      cvc-mesher/LBIE/LBIE_geoframe.cpp
      cvc-mesher/LBIE/tetra.cpp
      cvc-mesher/LBIE/hexa.cpp
      cvc-mesher/LBIE/LBIE_Mesher.cpp
      #cvc-mesher/LBIE/test.cpp
      cvc-mesher/LBIE/octree.cpp
      #cvc-mesher/Mesher/main.cpp
  )
endif(CVC_ENABLE_MESHER)

add_library(cvc SHARED ${SOURCE_FILES} ${INCLUDE_FILES})

SetupBoost(cvc)

target_link_libraries(cvc ${CVC_LINK_LIBS})

add_executable(trans-cvc trans-cvc.cpp)
target_link_libraries(trans-cvc cvc)

