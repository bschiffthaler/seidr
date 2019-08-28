# - Find CLP (includes and library)
# This module defines
#  CLP_INCLUDE_DIR
#  CLP_LIBRARIES
#  CLP_FOUND
# also defined, but not for general use are
#  CLP_LIBRARY, where to find the library.

if (NOT "$ENV{CLP_ROOT}" STREQUAL "")
  SET(CLP_CUSTOM_PATH $ENV{CLP_ROOT})
endif()

FIND_PATH(CLP_INCLUDE_DIR coin
/usr/include/
/usr/local/include/
/include
${CLP_CUSTOM_PATH}/include
)

SET(CLP_NAMES "Clp")
FIND_LIBRARY(CLP_LIBRARY
  NAMES ${CLP_NAMES}
  PATHS /usr/lib64 /usr/lib /usr/local/lib /lib /usr/lib/coin /usr/local/lib/coin /lib/coin ${CLP_CUSTOM_PATH}/lib
  )

SET(CLPSOLVER_NAMES "ClpSolver")
FIND_LIBRARY(CLPSOLVER_LIBRARY
  NAMES ${CLPSOLVER_NAMES}
  PATHS /usr/lib64 /usr/lib /usr/local/lib /lib /usr/lib/coin /usr/local/lib/coin /lib/coin ${CLP_CUSTOM_PATH}/lib
  )

SET(COINUTILS_NAMES "CoinUtils")
FIND_LIBRARY(COINUTILS_LIBRARY
  NAMES ${COINUTILS_NAMES}
  PATHS /usr/lib64 /usr/lib /usr/local/lib /lib /usr/lib/coin /usr/local/lib/coin /lib/coin ${CLP_CUSTOM_PATH}/lib
  )

IF (CLP_LIBRARY AND CLPSOLVER_LIBRARY AND COINUTILS_LIBRARY AND CLP_INCLUDE_DIR)
    SET(CLP_LIBRARIES ${CLP_LIBRARY} ${CLPSOLVER_LIBRARY} ${COINUTILS_LIBRARY})
    SET(CLP_FOUND "YES")
ELSE ()
  SET(CLP_FOUND "NO")
ENDIF ()


IF (CLP_FOUND)
   IF (NOT CLP_FIND_QUIETLY)
      MESSAGE(STATUS "Found CLP libraries: ${CLP_LIBRARIES}")
   ENDIF (NOT CLP_FIND_QUIETLY)
ELSE (CLP_FOUND)
   IF (CLP_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find a CLP library")
   ENDIF (CLP_FIND_REQUIRED)
ENDIF (CLP_FOUND)

# Deprecated declarations.
SET (NATIVE_CLP_INCLUDE_PATH ${CLP_INCLUDE_DIR} )
GET_FILENAME_COMPONENT (NATIVE_CLP_LIB_PATH ${CLP_LIBRARY} PATH)

MARK_AS_ADVANCED(
  CLP_LIBRARY
  CLP_INCLUDE_DIR
  )
