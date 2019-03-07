# - Find TBB (includes and library)
# This module defines
#  TBB_INCLUDE_DIR
#  TBB_LIBRARIES
#  TBB_FOUND
# also defined, but not for general use are
#  TBB_LIBRARY, where to find the library.

if (NOT "$ENV{TBB_ROOT}" STREQUAL "")
  SET(TBB_CUSTOM_PATH $ENV{TBB_ROOT})
endif()

FIND_PATH(TBB_INCLUDE_DIR tbb/task.h
/usr/include/
/usr/local/include/
/include
${TBB_CUSTOM_PATH}/include
)

SET(TBB_NAMES ${TBB_NAMES} tbb)
FIND_LIBRARY(TBB_LIBRARY
  NAMES ${TBB_NAMES}
  PATHS /usr/lib /usr/local/lib /lib /usr/lib/tbb /usr/local/lib/tbb /lib/tbb ${TBB_CUSTOM_PATH}/lib
  )

IF (TBB_LIBRARY AND TBB_INCLUDE_DIR)
    SET(TBB_LIBRARIES ${TBB_LIBRARY})
    SET(TBB_FOUND "YES")
ELSE (TBB_LIBRARY AND TBB_INCLUDE_DIR)
  SET(TBB_FOUND "NO")
ENDIF (TBB_LIBRARY AND TBB_INCLUDE_DIR)


IF (TBB_FOUND)
   IF (NOT TBB_FIND_QUIETLY)
      MESSAGE(STATUS "Found a TBB library: ${TBB_LIBRARIES}")
   ENDIF (NOT TBB_FIND_QUIETLY)
ELSE (TBB_FOUND)
   IF (TBB_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find a TBB library")
   ENDIF (TBB_FIND_REQUIRED)
ENDIF (TBB_FOUND)

# Deprecated declarations.
SET (NATIVE_TBB_INCLUDE_PATH ${TBB_INCLUDE_DIR} )
GET_FILENAME_COMPONENT (NATIVE_TBB_LIB_PATH ${TBB_LIBRARY} PATH)

MARK_AS_ADVANCED(
  TBB_LIBRARY
  TBB_INCLUDE_DIR
  )
