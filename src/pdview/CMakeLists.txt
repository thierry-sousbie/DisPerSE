set(C_DIR "${PROJECT_SOURCE_DIR}/src/C/")
set(MSE_DIR "${PROJECT_SOURCE_DIR}/src/mse/")
include_directories(${C_DIR})
include_directories(${MSE_DIR})

if(MATHGL_FOUND AND Qt5_FOUND)
  
  SET( PDVIEW_SRCS
    pdview.cpp
    ${MSE_DIR}/NDskel_breakdown.cpp 
    ${C_DIR}/NDnetwork.c
    ${C_DIR}/NDnetwork_tags.c
    ${C_DIR}/NDnet_VTK_IO.c
    ${C_DIR}/NDnet_PLY_IO.c 
    ${C_DIR}/rply/rply.c
    ${C_DIR}/gadget_io.c 
    ${C_DIR}/NDskeleton.c 
    ${C_DIR}/NDskel_tags.c 
    ${C_DIR}/NDskel_VTK_IO.c
    ${C_DIR}/NDfield.c 
    ${C_DIR}/NDfield_VTK_IO.c 
    ${C_DIR}/mystring.c
    ${C_DIR}/myendian.c
    #${C_DIR}/sort.c
    ${MSE_DIR}/sortID.cpp
    ${C_DIR}/simplex.c
    ${C_DIR}/probas.c
    )
  
  SET( PDVIEW_MOC_HDRS
    PDV_interface.hxx
    PD_widget.hxx
    )
  
  SET( PDVIEW_UIS
    
    )
  
  SET( PDVIEW_RCS
    help.qrc
    )
  
  # enable warnings
  ADD_DEFINITIONS( -Wall )
  qt5_add_resources(PDVIEW_RC_SRCS ${PDVIEW_RCS})
  qt5_wrap_ui(PDVIEW_UI_HDRS ${PDVIEW_UIS})
  qt5_wrap_cpp(PDVIEW_MOC_SRCS ${PDVIEW_MOC_HDRS})
  INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR})
  
  ADD_EXECUTABLE( pdview ${PDVIEW_SRCS} ${PDVIEW_MOC_SRCS} ${PDVIEW_RC_SRCS} ${PDVIEW_UI_HDRS} )
    # Add the paths to certain libraries in case we are in a conda env
    set(MATHGL_LIBRARIES 
	 "$ENV{CONDA_PREFIX}/lib/libmgl-qt5.so"
         "$ENV{CONDA_PREFIX}/lib/libmgl.so")
    # Link Qt5 and MathGL
    target_link_libraries(pdview
        Qt5::Core
        Qt5::Widgets
        Qt5::Gui
        Qt5::OpenGL
        ${MATHGL_LIBRARIES}
    )

  if(GSL_FOUND)
    target_link_libraries(pdview ${GSL_LIBRARIES})
  endif(GSL_FOUND)

  SET( TARGETS_LIST ${TARGETS_LIST} pdview PARENT_SCOPE)
  install (TARGETS pdview DESTINATION bin)
else()
  SET( NOT_BUILT_TARGETS_LIST ${NOT_BUILT_TARGETS_LIST} pdview PARENT_SCOPE)
  SET( DISABLED_FEATURES_LIST_LOC ${DISABLED_FEATURES_LIST_LOC} mse:interactive)

endif(MATHGL_QT_FOUND AND QT_FOUND)

SET(DISABLED_FEATURES_LIST ${DISABLED_FEATURES_LIST} ${DISABLED_FEATURES_LIST_LOC} PARENT_SCOPE)
