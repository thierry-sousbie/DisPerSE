# Try to find gnu scientific library SDL
# See 
# http://www.gnu.org/software/sdl/  and 
# http://gnuwin32.sourceforge.net/packages/sdl.htm
#
# Once run this will define: 
# 
# SDL_FOUND       = system has SDL lib
#
# SDL_LIBRARIES   = full path to the libraries
#    on Unix/Linux with additional linker flags from "sdl-config --libs"
# 
# CMAKE_SDL_CXX_FLAGS  = Unix compiler flags for SDL, essentially "`sdl-config --cxxflags`"
#
# SDL_INCLUDE_DIR      = where to find headers 
#
# SDL_LINK_DIRECTORIES = link directories, useful for rpath on Unix
# SDL_EXE_LINKER_FLAGS = rpath on Unix
#
# Felix Woelk 07/2004
# Jan Woetzel
#
# www.mip.informatik.uni-kiel.de
# --------------------------------

FIND_PROGRAM(SDL_CONFIG sdl-config
  ${SDL_DIR}/bin
  NO_DEFAULT_PATH
  )

SET(SDL_CONFIG_PREFER_PATH 
  "$ENV{SDL_DIR}/bin"
  "$ENV{SDL_DIR}"
  "$ENV{SDL_HOME}/bin" 
  "$ENV{SDL_HOME}" 
  "${SDL_DIR}/bin"
  CACHE STRING "preferred path to SDL (sdl-config)"
  )
FIND_PROGRAM(SDL_CONFIG sdl-config
  ${SDL_CONFIG_PREFER_PATH}
  /usr/bin/
  )


IF (SDL_CONFIG) 
  MESSAGE(STATUS "SDL using sdl-config ${SDL_CONFIG}")
      # set CXXFLAGS to be fed into CXX_FLAGS by the user:
      EXEC_PROGRAM(${SDL_CONFIG}
        ARGS --cflags
        OUTPUT_VARIABLE  SDL_CXX_FLAGS )
      #SET(SDL_CXX_FLAGS "`${SDL_CONFIG} --cflags`")
      
      # set INCLUDE_DIRS to prefix+include
      EXEC_PROGRAM(${SDL_CONFIG}
        ARGS --prefix
        OUTPUT_VARIABLE SDL_PREFIX)
      SET(SDL_INCLUDE_DIR ${SDL_PREFIX}/include CACHE STRING INTERNAL)

      # set link libraries and link flags
      
      #SET(SDL_LIBRARIES "`${SDL_CONFIG} --libs`")
      
      # extract link dirs for rpath  
      EXEC_PROGRAM(${SDL_CONFIG}
        ARGS --libs
        OUTPUT_VARIABLE  SDL_CONFIG_LIBS )
      SET(SDL_LIBRARIES "${SDL_CONFIG_LIBS} -lSDL_image")

      # split off the link dirs (for rpath)
      # use regular expression to match wildcard equivalent "-L*<endchar>"
      # with <endchar> is a space or a semicolon
      STRING(REGEX MATCHALL "[-][L]([^ ;])+" 
        SDL_LINK_DIRECTORIES_WITH_PREFIX 
        "${SDL_CONFIG_LIBS}" )
      #      MESSAGE("DBG  SDL_LINK_DIRECTORIES_WITH_PREFIX=${SDL_LINK_DIRECTORIES_WITH_PREFIX}")

      # remove prefix -L because we need the pure directory for LINK_DIRECTORIES
      
      IF (SDL_LINK_DIRECTORIES_WITH_PREFIX)
        STRING(REGEX REPLACE "[-][L]" "" SDL_LINK_DIRECTORIES ${SDL_LINK_DIRECTORIES_WITH_PREFIX} )
      ENDIF (SDL_LINK_DIRECTORIES_WITH_PREFIX)
      SET(SDL_EXE_LINKER_FLAGS "-Wl,-rpath,${SDL_LINK_DIRECTORIES}" CACHE STRING INTERNAL)
      #      MESSAGE("DBG  SDL_LINK_DIRECTORIES=${SDL_LINK_DIRECTORIES}")
      #      MESSAGE("DBG  SDL_EXE_LINKER_FLAGS=${SDL_EXE_LINKER_FLAGS}")

      #      ADD_DEFINITIONS("-DHAVE_SDL")
      #      SET(SDL_DEFINITIONS "-DHAVE_SDL")
      MARK_AS_ADVANCED(
        SDL_CXX_FLAGS
        SDL_INCLUDE_DIR
        SDL_LIBRARIES
        SDL_LINK_DIRECTORIES
        SDL_DEFINITIONS
	)
      MESSAGE(STATUS "Using SDL from ${SDL_PREFIX}")
ELSE(SDL_CONFIG)
  INCLUDE(UsePkgConfig) #needed for PKGCONFIG(...)
  
  MESSAGE(STATUS "SDL using pkgconfig")
  PKGCONFIG(sdl SDL_INCLUDE_DIR SDL_LINK_DIRECTORIES SDL_LIBRARIES SDL_CXX_FLAGS)
  IF(SDL_INCLUDE_DIR)
    MARK_AS_ADVANCED(
      SDL_CXX_FLAGS
      SDL_INCLUDE_DIR
      SDL_LIBRARIES
      SDL_LINK_DIRECTORIES
      )
    
  ELSE(SDL_INCLUDE_DIR)      
    MESSAGE("FindSDL.cmake: sdl-config/pkg-config sdl not found. Please set it manually. SDL_DIR=path_to_sdl/")
  ENDIF(SDL_INCLUDE_DIR)
  
ENDIF(SDL_CONFIG)

IF(SDL_LIBRARIES)
  IF(SDL_INCLUDE_DIR OR SDL_CXX_FLAGS)

    SET(SDL_FOUND 1)
    IF(EXISTS "${SDL_INCLUDE_DIR}/SDL/SDL_image.h")
      SET(SDL-IMAGE_FOUND 1)
    ENDIF()
    
  ENDIF(SDL_INCLUDE_DIR OR SDL_CXX_FLAGS)
ENDIF(SDL_LIBRARIES)


# ==========================================
IF(NOT SDL_FOUND)
  # make FIND_PACKAGE friendly
  IF(NOT SDL_FIND_QUIETLY)
    IF(SDL_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "SDL required, please specify it's location.")
    ELSE(SDL_FIND_REQUIRED)
      MESSAGE(STATUS       "WARNING: SDL was not found.")
    ENDIF(SDL_FIND_REQUIRED)
  ENDIF(NOT SDL_FIND_QUIETLY)
ELSEIF (NOT SDL-IMAGE_FOUND)
  IF(NOT SDL_FIND_QUIETLY)
    IF(SDL_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "ERROR: SDL-image required, please specify it's location. (SDL was found)")
    ELSE(SDL_FIND_REQUIRED)
      MESSAGE(STATUS       "WARNING: SDL was found but not SDL-image.")
    ENDIF(SDL_FIND_REQUIRED)
  ENDIF(NOT SDL_FIND_QUIETLY)
ENDIF(NOT SDL_FOUND)
