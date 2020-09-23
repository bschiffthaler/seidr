cmake_minimum_required(VERSION 3.9.4)
project(deps)

if(CMAKE_BUILD_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(ExternalProject)

ExternalProject_Add(
  NetworKit
  CMAKE_ARGS
    "-DNETWORKIT_STATIC=ON"
    "-DCMAKE_CXX_FLAGS='-include stdexcept'"
    "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/NetworKit"
  GIT_REPOSITORY https://github.com/networkit/networkit.git
  GIT_SHALLOW 1
  GIT_PROGRESS 1
  GIT_TAG ${NETWORKIT_VERSION}
  UPDATE_DISCONNECTED ON
)

# Build ranger as a library
file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/include/ranger)
file(GLOB RANGER_HEADERS "../external/ranger/src/*.h")
file(COPY ${RANGER_HEADERS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/ranger)
add_library(ranger ../external/ranger/src/Data.cpp
  ../external/ranger/src/Data.cpp
  ../external/ranger/src/Forest.cpp ../external/ranger/src/ForestRegression.cpp 
  ../external/ranger/src/Tree.cpp ../external/ranger/src/TreeRegression.cpp 
  ../external/ranger/src/utility.cpp)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(ranger PRIVATE ${ARMADILLO_LIBRARIES})

# Build libbs
file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/include/libbs)
file(GLOB BS_HEADERS "../external/libbs/src/*.h")
file(COPY ${BS_HEADERS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/bs)
add_library(bs ../external/libbs/src/vitter_a.cpp ../external/libbs/src/vitter_d.cpp
  ../external/libbs/src/aux.cpp ../external/libbs/src/describe.cpp 
  ../external/libbs/src/histogram.cpp ../external/libbs/src/common.cpp) 


