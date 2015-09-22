# - FindMathGL.cmake
# This module can be used to find MathGL and several of its optional components.
#
# You can specify one or more component as you call this find module.
# Possible components are: FLTK, GLUT, Qt, WX.
#
# The following variables will be defined for your use:
#
#  MATHGL_FOUND           = MathGL and all specified components found
#  MATHGL_INCLUDE_DIRS    = The MathGL include directories
#  MATHGL_LIBRARIES       = The libraries to link against to use MathGL
#                           and all specified components
#  MATHGL_VERSION_STRING  = A human-readable version of the MathGL (e.g. 1.11)
#  MATHGL_XXX_FOUND       = Component XXX found (replace XXX with uppercased
#                           component name -- for example, QT or FLTK)
#
# The minimum required version and needed components can be specified using
# the standard find_package()-syntax, here are some examples:
#  find_package(MathGL 1.11 Qt REQUIRED) - 1.11 + Qt interface, required
#  find_package(MathGL 1.10 REQUIRED)    - 1.10 (no interfaces), required
#  find_package(MathGL 1.10 Qt WX)       - 1.10 + Qt and WX interfaces, optional
#  find_package(MathGL 1.11)             - 1.11 (no interfaces), optional
#
# Typical usage could be something like this:
#   find_package(MathGL 1.11 GLUT REQUIRED)
#   include_directories(${MATHGL_INCLUDE_DIRS})
#   add_executable(myexe main.cpp)
#   target_link_libraries(myexe ${MATHGL_LIBRARIES})
#

#=============================================================================
# Copyright (c) 2011 Denis Pesotsky <denis@kde.ru>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file COPYING-CMAKE-MODULES for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

FIND_PATH(MATHGL_INCLUDE_DIR
  NAMES mgl/mgl.h 
  PATHS ${MATHGL_DIR}/include/
  NO_DEFAULT_PATH
  )

FIND_PATH(MATHGL_INCLUDE_DIR
          NAMES mgl/mgl.h 
	  PATHS ${MATHGL_DIR}/include/
          DOC "The MathGL include directory")

FIND_LIBRARY(MATHGL_LIBRARY
             NAMES mgl
             PATHS ${MATHGL_DIR}/include/
             NO_DEFAULT_PATH)

FIND_LIBRARY(MATHGL_LIBRARY
             NAMES mgl
             PATHS ${MATHGL_DIR}/lib/ ${MATHGL_LIBRARY_DIR}
             DOC "The MathGL LIBRARY directory")

GET_FILENAME_COMPONENT(MATHGL_LIBRARY_DIR ${MATHGL_LIBRARY} PATH)

SET(MATHGL_LIBRARIES ${MATHGL_LIBRARY})
SET(MATHGL_INCLUDE_DIRS ${MATHGL_INCLUDE_DIR})

IF(MATHGL_INCLUDE_DIR)
  SET(_CONFIG_FILE_NAME "mgl/config.h")
  SET(_CONFIG_FILE_PATH "${MATHGL_INCLUDE_DIR}/${_CONFIG_FILE_NAME}")
  SET(_VERSION_ERR "Cannot determine MathGL version")
  IF(EXISTS "${_CONFIG_FILE_PATH}")
    FILE(STRINGS "${_CONFIG_FILE_PATH}"
         MATHGL_VERSION_STRING REGEX "^#define PACKAGE_VERSION \"[^\"]*\"$")
    IF(MATHGL_VERSION_STRING)
      STRING(REGEX
             REPLACE "^#define PACKAGE_VERSION \"([^\"]*)\"$" "\\1"
             MATHGL_VERSION_STRING ${MATHGL_VERSION_STRING})
    ELSE()
      MESSAGE(FATAL_ERROR "${_VERSION_ERR}: ${_CONFIG_FILE_NAME} parse error")
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR "${_VERSION_ERR}: ${_CONFIG_FILE_NAME} not found")
  ENDIF()
ENDIF()

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MathGL
                                  REQUIRED_VARS MATHGL_LIBRARY
                                                MATHGL_INCLUDE_DIR
                                  VERSION_VAR MATHGL_VERSION_STRING)

FOREACH(_Component ${MathGL_FIND_COMPONENTS})
  STRING(TOLOWER ${_Component} _component)
  STRING(TOUPPER ${_Component} _COMPONENT)
  
  SET(MathGL_${_Component}_FIND_REQUIRED ${MathGL_FIND_REQUIRED})
  SET(MathGL_${_Component}_FIND_QUIETLY true)
  
  FIND_PATH(MATHGL_${_COMPONENT}_INCLUDE_DIR
            NAMES mgl/mgl_${_component}.h
            PATHS ${MATHGL_INCLUDE_DIR} NO_DEFAULT_PATH)
  FIND_LIBRARY(MATHGL_${_COMPONENT}_LIBRARY
               NAMES mgl-${_component}
               PATHS ${MATHGL_LIBRARY_DIR} NO_DEFAULT_PATH)

  FIND_PACKAGE_HANDLE_STANDARD_ARGS(MathGL_${_Component} DEFAULT_MSG
                                    MATHGL_${_COMPONENT}_LIBRARY
                                    MATHGL_${_COMPONENT}_INCLUDE_DIR)
  
  IF(MATHGL_${_COMPONENT}_FOUND)
    SET(MATHGL_LIBRARIES
        ${MATHGL_LIBRARIES} ${MATHGL_${_COMPONENT}_LIBRARY})
    SET(MATHGL_INCLUDE_DIRS
        ${MATHGL_INCLUDE_DIRS} ${MATHGL_${_COMPONENT}_INCLUDE_DIR})
  ENDIF()

  MARK_AS_ADVANCED(MATHGL_${_COMPONENT}_INCLUDE_DIR
                   MATHGL_${_COMPONENT}_LIBRARY)
ENDFOREACH()

MARK_AS_ADVANCED(MATHGL_INCLUDE_DIR MATHGL_LIBRARY)
