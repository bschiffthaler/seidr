cmake_minimum_required(VERSION 3.9.4)
project(libnetworkit)

if(CMAKE_BUILD_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif()

# Build Networkit as library 
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/include/networkit)
file(COPY ../external/networkit/networkit/cpp DESTINATION 
     ${CMAKE_CURRENT_BINARY_DIR}/include)
file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/include/cpp
     ${CMAKE_CURRENT_BINARY_DIR}/include/networkit)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
add_library(networkit 
  ../external/networkit/networkit/cpp/base/Algorithm.cpp
  ../external/networkit/networkit/cpp/centrality/ApproxCloseness.cpp
  ../external/networkit/networkit/cpp/centrality/Betweenness.cpp
  ../external/networkit/networkit/cpp/centrality/Centrality.cpp
  ../external/networkit/networkit/cpp/centrality/Closeness.cpp
  ../external/networkit/networkit/cpp/centrality/EigenvectorCentrality.cpp
  ../external/networkit/networkit/cpp/centrality/EstimateBetweenness.cpp
  ../external/networkit/networkit/cpp/centrality/KatzCentrality.cpp
  ../external/networkit/networkit/cpp/centrality/LocalClusteringCoefficient.cpp
  ../external/networkit/networkit/cpp/centrality/PageRank.cpp
  ../external/networkit/networkit/cpp/centrality/SpanningEdgeCentrality.cpp
  ../external/networkit/networkit/cpp/components/ConnectedComponents.cpp
  ../external/networkit/networkit/cpp/components/ParallelConnectedComponents.cpp
  ../external/networkit/networkit/cpp/distance/Diameter.cpp
  ../external/networkit/networkit/cpp/global/ClusteringCoefficient.cpp
  ../external/networkit/networkit/cpp/graph/Graph.cpp
  ../external/networkit/networkit/cpp/graph/GraphBuilder.cpp
  ../external/networkit/networkit/cpp/auxiliary/Log.cpp
  ../external/networkit/networkit/cpp/auxiliary/NumericTools.cpp
  ../external/networkit/networkit/cpp/auxiliary/Random.cpp
  ../external/networkit/networkit/cpp/auxiliary/SignalHandling.cpp
  ../external/networkit/networkit/cpp/auxiliary/Timer.cpp
  ../external/networkit/networkit/cpp/algebraic/CSRMatrix.cpp
  ../external/networkit/networkit/cpp/algebraic/DenseMatrix.cpp
  ../external/networkit/networkit/cpp/algebraic/Vector.cpp
  ../external/networkit/networkit/cpp/coarsening/ParallelPartitionCoarsening.cpp
  ../external/networkit/networkit/cpp/distance/BFS.cpp
  ../external/networkit/networkit/cpp/distance/Dijkstra.cpp
  ../external/networkit/networkit/cpp/distance/SSSP.cpp
  ../external/networkit/networkit/cpp/distance/APSP.cpp
  ../external/networkit/networkit/cpp/coarsening/GraphCoarsening.cpp
  ../external/networkit/networkit/cpp/structures/Partition.cpp
  )

# Build ranger as a library
file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/include/ranger)
file(GLOB RANGER_HEADERS "../external/ranger/src/*.h")
file(COPY ${RANGER_HEADERS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include)
add_library(ranger ../external/ranger/src/Data.cpp ../external/ranger/src/DataDouble.cpp 
  ../external/ranger/src/DataChar.cpp ../external/ranger/src/DataFloat.cpp 
  ../external/ranger/src/Forest.cpp ../external/ranger/src/ForestRegression.cpp 
  ../external/ranger/src/Tree.cpp ../external/ranger/src/TreeRegression.cpp 
  ../external/ranger/src/utility.cpp)


