cmake_minimum_required(VERSION 3.9.4)
project(SeidrFile)

set(CMAKE_CXX_STANDARD 11)

include_directories(../../external/bgzf ../../src/Seidr)

add_library(seidrfile SHARED 
            ../../src/Seidr/Serialize.cpp
            ../../external/bgzf/bgzf.c
            src/wrapper.cpp)

set_target_properties(seidrfile PROPERTIES SUFFIX ".so")

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( seidrfile ${ZLIB_LIBRARIES} )

find_package( Boost COMPONENTS system filesystem python37 REQUIRED )
include_directories( ${BOOST_INCLUDE_DIR} )
target_link_libraries( seidrfile ${Boost_LIBRARIES})

find_package (Python COMPONENTS Development REQUIRED)
include_directories( ${Python_INCLUDE_DIRS} )
target_link_libraries( seidrfile ${Python_LIBRARIES})