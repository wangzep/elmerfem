# CMake script for finding Permon for Elmer
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(_permonIfaceF90 "permon_iface.f90" INTERNAL)
SET(_permonLibName "permon" INTERNAL)

# If MKL_LIBRARIES libraries are already defined, do nothing
SET(Permon_FOUND FALSE)

SET(_permonIncludePaths
  "$ENV{PERMON_ROOT}/include"
  "${PERMON_ROOT}/include"
  INTERNAL
)

SET(_permonLibPaths
  "$ENV{PERMON_ROOT}/lib"
  "${PERMON_ROOT}/lib"
  INTERNAL
)

SET(_permonInterfaceSrcPaths
  "ENV{PERMON_ROOT}/lib/${_pmonIfaceF90}"
  "${PERMON_ROOT}/lib/${_pmonIfaceF90}"
  INTERNAL
)

# Find Permon library
#FIND_LIBRARY(PERMON_LIBRARIES ${_permonLibName}${SHL_EXTENSION} HINTS ${_permonLibPaths})

SET(PERMON_LIBRARIES /home/silvonen/work/permon_test/libpermon.a CACHE FILE "")
SET(PERMON_INTERFACE_SOURCE /home/silvonen/work/permon_test/foo.F90 CACHE FILE "")

# Find the actual interface file
#FIND_FILE(PERMON_INTERFACE_SOURCE NAMES ${_permonIfaceF90} PATHS ${_permonInterfaceSrcPaths})

#IF(_permonLibrary AND _permonInterfaceSrc)
#  UNSET(PERMON_FAILMSG)
#  SET(PERMON_LIBRARIES ${_permonLibrary} CACHE FILEPATH "Permon library")
#  SET(PERMON_INTERFACE_SOURCE ${_permonInterfaceSrc} CACHE FILEPATH "Permon interface source")
#ELSE()
#  SET(PERMON_FAILMSG "Permon not found")
#ENDIF()

IF(PERMON_LIBRARIES AND PERMON_INTERFACE_SOURCE)
  SET(Permon_FOUND TRUE)
ENDIF()

IF(Permon_FOUND)
  IF (NOT PERMON_FIND_QUIETLY)
    MESSAGE(STATUS "A library with Permon API found.")
  ENDIF()
ELSE()
  IF (PERMON_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${PERMON_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  PERMON_LIBRARIES
  PERMON_INTERFACE_SOURCE
)
