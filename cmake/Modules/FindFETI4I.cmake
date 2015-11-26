# CMake script for finding Permon for Elmer
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(_feti4iIfaceF90 "feti4i_mod.F90")
SET(_feti4iLibName "libfeti4i.so")

# If MKL_LIBRARIES libraries are already defined, do nothing
SET(FETI4I_FOUND FALSE)

SET(_feti4iIncludePaths
  "$ENV{FETI4I_ROOT}/include"
  "${FETI4I_ROOT}/include"
  INTERNAL
)

SET(_feti4iLibPaths
  "$ENV{FETI4I_ROOT}/lib"
  "${FETI4I_ROOT}/lib"
  INTERNAL
)

SET(_feti4iInterfaceSrcPaths
  "ENV{FETI4I_ROOT}/lib/${_feti4iIfaceF90}"
  "${FETI4I_ROOT}/lib/${_feti4iIfaceF90}"
  INTERNAL
)

# Find Feti4i library
#FIND_LIBRARY(FETI4I_LIBRARIES ${_feti4iLibName}${SHL_EXTENSION} HINTS ${_feti4iLibPaths})

SET(FETI4I_LIBRARIES ${PROJECT_SOURCE_DIR}/feti4i/feti4i_fortran_test/lib/${_feti4iLibName} CACHE FILE "")
SET(FETI4I_INTERFACE_SOURCE ${PROJECT_SOURCE_DIR}/feti4i/${_feti4iIfaceF90} CACHE FILE "")

# Find the actual interface file
#FIND_FILE(FETI4I_INTERFACE_SOURCE NAMES ${_feti4iIfaceF90} PATHS ${_feti4iInterfaceSrcPaths})

#IF(_feti4iLibrary AND _feti4iInterfaceSrc)
#  UNSET(FETI4I_FAILMSG)
#  SET(FETI4I_LIBRARIES ${_feti4iLibrary} CACHE FILEPATH "Feti4i library")
#  SET(FETI4I_INTERFACE_SOURCE ${_feti4iInterfaceSrc} CACHE FILEPATH "Feti4i interface source")
#ELSE()
#  SET(FETI4I_FAILMSG "Feti4i not found")
#ENDIF()

IF(FETI4I_LIBRARIES AND FETI4I_INTERFACE_SOURCE)
  SET(FETI4I_FOUND TRUE)
ENDIF()

IF(FETI4I_FOUND)
  IF (NOT FETI4I_FIND_QUIETLY)
    MESSAGE(STATUS "A library with FETI4I API found.")
  ENDIF()
ELSE()
  IF (FETI4I_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${FETI4I_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  FETI4I_LIBRARIES
  FETI4I_INTERFACE_SOURCE
)
