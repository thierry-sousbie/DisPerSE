set(C_DIR "${PROJECT_SOURCE_DIR}/src/C/")
set(MSE_DIR "${PROJECT_SOURCE_DIR}/src/mse/")
include_directories(${C_DIR})
include_directories(${MSE_DIR})

if(GSL_FOUND AND CGAL_FOUND)

  include( ${CGAL_USE_FILE} )

  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    STRING(REGEX REPLACE "-O2" " "  CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
  endif()

  #message(STATUS "  USED   ${CMAKE_CXX_FLAGS}")

  SET(DELAUNAY2D_SRC
    delaunay_2D.cpp 
    ${MSE_DIR}/NDskel_breakdown.cpp 
#    ${C_DIR}/pthreadWrapper.c    
    ${C_DIR}/mystring.c
    ${C_DIR}/gadget_io.c 
    ${C_DIR}/NDfield.c 
    ${C_DIR}/NDfield_VTK_IO.c 
    ${C_DIR}/myendian.c 
    ${C_DIR}/simplex.c 
    ${C_DIR}/NDnetwork.c    
    ${C_DIR}/NDnetwork_tags.c 
    ${C_DIR}/NDnet_VTK_IO.c
    ${C_DIR}/NDnet_PLY_IO.c 
    ${C_DIR}/NDskeleton.c 
    ${C_DIR}/NDskel_tags.c 
    ${C_DIR}/NDskel_VTK_IO.c  
    ${C_DIR}/rply/rply.c
    ${C_DIR}/dtfi.c 
    ${C_DIR}/dtfi.h
    #${C_DIR}/sort.c 
    ${C_DIR}/distances.c
    ${C_DIR}/asciiSurvey.c 
    ${C_DIR}/healpix.c
    ${MSE_DIR}/sortID.cpp
    )

  SET( DELAUNAY3D_SRC
    delaunay_3D.cpp 
    ${MSE_DIR}/NDskel_breakdown.cpp 
#    ${C_DIR}/pthreadWrapper.c    
    ${C_DIR}/mystring.c
    ${C_DIR}/gadget_io.c 
    ${C_DIR}/NDfield.c
    ${C_DIR}/NDfield_VTK_IO.c 
    ${C_DIR}/myendian.c 
    ${C_DIR}/simplex.c 
    ${C_DIR}/NDnetwork.c    
    ${C_DIR}/NDnetwork_tags.c 
    ${C_DIR}/NDnet_VTK_IO.c 
    ${C_DIR}/NDnet_PLY_IO.c 
    ${C_DIR}/NDskeleton.c 
    ${C_DIR}/NDskel_tags.c 
    ${C_DIR}/NDskel_VTK_IO.c   
    ${C_DIR}/rply/rply.c
    ${C_DIR}/dtfi.c 
    ${C_DIR}/dtfi.h
    #${C_DIR}/sort.c
    ${C_DIR}/distances.c
    ${C_DIR}/asciiSurvey.c 
    ${C_DIR}/healpix.c
    ${MSE_DIR}/sortID.cpp
    )

  #include( CGAL_CreateSingleSourceCGALProgram )
  #create_single_source_cgal_program(${DELAUNAY2D_SRC})
  #create_single_source_cgal_program(${DELAUNAY3D_SRC})

  add_executable(delaunay_2D ${DELAUNAY2D_SRC})
  add_executable(delaunay_3D ${DELAUNAY3D_SRC})
  
  #target_link_libraries(delaunay_2D ${GSL_LIBRARIES})
  #target_link_libraries(delaunay_3D ${GSL_LIBRARIES})

  if (CFITSIO_FOUND)
    #target_link_libraries(delaunay_2D ${CFITSIO_LIBRARIES})
    #target_link_libraries(delaunay_3D ${CFITSIO_LIBRARIES})
  else()
    SET( DISABLED_FEATURES_LIST_LOC ${DISABLED_FEATURES_LIST_LOC} delaunay::angular_mask)
  endif()

  if (SDL-IMAGE_FOUND)
    #target_link_libraries(delaunay_2D ${SDL_LIBRARIES})
    #target_link_libraries(delaunay_3D ${SDL_LIBRARIES})
  endif()

  SET( TARGETS_LIST ${TARGETS_LIST} delaunay_2D delaunay_3D PARENT_SCOPE)
  install (TARGETS delaunay_2D delaunay_3D DESTINATION bin)
else()
  SET( NOT_BUILT_TARGETS_LIST ${NOT_BUILT_TARGETS_LIST} delaunay_2D delaunay_3D PARENT_SCOPE)
endif(GSL_FOUND AND CGAL_FOUND)

SET(DISABLED_FEATURES_LIST ${DISABLED_FEATURES_LIST} ${DISABLED_FEATURES_LIST_LOC} PARENT_SCOPE)

  
