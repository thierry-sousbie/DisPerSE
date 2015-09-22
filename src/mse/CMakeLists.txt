set(C_DIR "${PROJECT_SOURCE_DIR}/src/C/")
include_directories(${C_DIR})

if(GSL_FOUND)
  SET( MSE_SRC
    mse.cpp
    arcs_nodes.cpp
    NDskel_breakdown.cpp 
#    ${C_DIR}/pthreadWrapper.c    
    ${C_DIR}/NDnetwork.c     
    ${C_DIR}/NDnetwork_tags.c
    ${C_DIR}/NDnet_VTK_IO.c     
    ${C_DIR}/NDnet_PLY_IO.c 
    ${C_DIR}/rply/rply.c
    ${C_DIR}/NDskeleton.c
    ${C_DIR}/NDskel_tags.c
    ${C_DIR}/NDskel_VTK_IO.c 
    ${C_DIR}/NDfield.c
    ${C_DIR}/NDfield_VTK_IO.c 
    ${C_DIR}/distances.c
    ${C_DIR}/asciiSurvey.c
    ${C_DIR}/healpix.c
    ${C_DIR}/gadget_io.c 
    ${C_DIR}/fof_io.c 
    ${C_DIR}/simplex.c 
    ${C_DIR}/hypercube.c   
    #${C_DIR}/sort.c
    sortID.cpp
    ${C_DIR}/probas.c
    ${C_DIR}/unionfind.c 
    ${C_DIR}/smooth.c
    ${C_DIR}/mystring.c        
    ${C_DIR}/myendian.c 
    )
    
  add_executable(mse ${MSE_SRC})
  #target_link_libraries(mse ${GSL_LIBRARIES})
  
  if (SDL-IMAGE_FOUND)
    #target_link_libraries(mse ${SDL_LIBRARIES})
  else()
    SET(DISABLED_FEATURES_LIST_LOC ${DISABLED_FEATURES_LIST_LOC} mse:jpeg)
  endif()

  if (CFITSIO_FOUND)
    #target_link_libraries(mse ${CFITSIO_LIBRARIES})
  else()
    SET(DISABLED_FEATURES_LIST_LOC ${DISABLED_FEATURES_LIST_LOC} mse:FitsIO mse:Healpix)
  endif()

  SET(TARGETS_LIST_LOC ${TARGETS_LIST_LOC} mse)
  install (TARGETS mse DESTINATION bin)
else()
  SET(NOT_BUILT_TARGETS_LIST_LOC ${NOT_BUILT_TARGETS_LIST_LOC} mse)
endif(GSL_FOUND)

SET( SKELCONV_SRC
  skelconv.cpp 
  NDskel_breakdown.cpp 
#  ${C_DIR}/pthreadWrapper.c    
  ${C_DIR}/mystring.c 
  ${C_DIR}/NDskeleton.c 
  ${C_DIR}/NDskel_tags.c 
  ${C_DIR}/NDskel_VTK_IO.c 
  ${C_DIR}/NDnetwork.c
  ${C_DIR}/NDnetwork_tags.c  
  #${C_DIR}/sort.c     
  sortID.cpp
  ${C_DIR}/simplex.c 
  ${C_DIR}/smooth.c 
  ${C_DIR}/NDfield.c 
  ${C_DIR}/NDfield_VTK_IO.c 
  ${C_DIR}/myendian.c  
  ${C_DIR}/distances.c
  ${C_DIR}/asciiSurvey.c
  )

add_executable(skelconv ${SKELCONV_SRC})
if (CFITSIO_FOUND)
  #target_link_libraries(skelconv ${CFITSIO_LIBRARIES})
endif()

#if(GSL_FOUND)
#  target_link_libraries(skelconv ${GSL_LIBRARIES})
#endif(GSL_FOUND)

SET( TARGETS_LIST_LOC ${TARGETS_LIST_LOC} skelconv)
install (TARGETS skelconv DESTINATION bin)

SET( NETCONV_SRC
  netconv.cpp 
  NDskel_breakdown.cpp 
#  ${C_DIR}/pthreadWrapper.c    
  ${C_DIR}/mystring.c 
  ${C_DIR}/NDskeleton.c 
  ${C_DIR}/NDskel_tags.c 
  ${C_DIR}/NDskel_VTK_IO.c 
  ${C_DIR}/NDnetwork.c 
  ${C_DIR}/NDnetwork_tags.c 
  ${C_DIR}/NDnet_VTK_IO.c
  ${C_DIR}/NDnet_PLY_IO.c 
  ${C_DIR}/rply/rply.c
  ${C_DIR}/gadget_io.c 
  ${C_DIR}/distances.c
  ${C_DIR}/asciiSurvey.c
  ${C_DIR}/smooth.c 
  ${C_DIR}/NDfield.c 
  ${C_DIR}/NDfield_VTK_IO.c 
  ${C_DIR}/myendian.c 
  #${C_DIR}/sort.c 
  sortID.cpp
  ${C_DIR}/simplex.c
  )

add_executable(netconv ${NETCONV_SRC})
if (CFITSIO_FOUND)
  #target_link_libraries(netconv ${CFITSIO_LIBRARIES})
endif()

#if(GSL_FOUND)
#  target_link_libraries(netconv ${GSL_LIBRARIES})
#endif(GSL_FOUND)

SET(TARGETS_LIST_LOC ${TARGETS_LIST_LOC} netconv)
install (TARGETS netconv DESTINATION bin)

SET( FIELDCONV_SRC
  fieldconv.cpp 
  ${C_DIR}/mystring.c 
  #${C_DIR}/NDskeleton.c 
  #${C_DIR}/NDskel_tags.c 
  ${C_DIR}/NDnetwork.c 
  #${C_DIR}/NDnetwork_tags.c 
  #${C_DIR}/NDnet_VTK_IO.c 
  #${C_DIR}/smooth.c 
#  ${C_DIR}/pthreadWrapper.c    
  ${C_DIR}/NDfield.c 
  ${C_DIR}/NDfield_VTK_IO.c 
  ${C_DIR}/myendian.c 
  #${C_DIR}/sort.c 
  sortID.cpp
  ${C_DIR}/distances.c
  ${C_DIR}/asciiSurvey.c
  ${C_DIR}/simplex.c
  )

add_executable(fieldconv ${FIELDCONV_SRC})
if (CFITSIO_FOUND)
  #target_link_libraries(fieldconv ${CFITSIO_LIBRARIES})
endif()

if (SDL-IMAGE_FOUND)
  #target_link_libraries(fieldconv ${SDL_LIBRARIES})
endif()

#if(GSL_FOUND)
#  target_link_libraries(fieldconv ${GSL_LIBRARIES})
#endif(GSL_FOUND)

SET(TARGETS_LIST_LOC ${TARGETS_LIST_LOC} fieldconv)
install (TARGETS fieldconv DESTINATION bin)

SET(TARGETS_LIST ${TARGETS_LIST} ${TARGETS_LIST_LOC} PARENT_SCOPE)
SET(NOT_BUILT_TARGETS_LIST ${NOT_BUILT_TARGETS_LIST} ${NOT_BUILT_TARGETS_LIST_LOC} PARENT_SCOPE)

SET(DISABLED_FEATURES_LIST ${DISABLED_FEATURES_LIST} ${DISABLED_FEATURES_LIST_LOC} PARENT_SCOPE)
