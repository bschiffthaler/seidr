cmake_minimum_required(VERSION 3.9.4)
project(seidr)

enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

include_directories(src/External src/Common src/Narromi src/Tigress src/SVM 
  src/ElNet src/Seidr src/Correlation src/MI src/PLSNET src/PCor 
  external/ranger/src src/GENIE3 src/Anoverence)

set(SEIDR_MAJOR_VERSION 0)
set(SEIDR_MINOR_VERSION 9)
set(SEIDR_PATCH_VERSION 0)
set(SEIDR_VERSION
  ${SEIDR_MAJOR_VERSION}.${SEIDR_MINOR_VERSION}.${SEIDR_PATCH_VERSION})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEIDR_VERSION=${SEIDR_VERSION} ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=${SEIDR_VERSION} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_32BIT_WORD -DARMA_DONT_PRINT_ERRORS")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -DDEBUG=1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
set(CMAKE_Fortran_FLAGS "-m64 -fdefault-real-8 -fdefault-double-8 -ffixed-form -O3")

if(CMAKE_BUILD_STATIC)
  set(BUILD_SHARED_LIBRARIES OFF)
  set(Boost_USE_STATIC_LIBS ON)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

if(CMAKE_BUILD_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif()

if(NOT SEIDR_SCORE_FLOAT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEIDR_SCORE_DOUBLE")
endif()

# Build Networkit as library 
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/include/networkit)
file(COPY external/networkit/networkit/cpp DESTINATION 
     ${CMAKE_CURRENT_BINARY_DIR}/include)
file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/include/cpp
     ${CMAKE_CURRENT_BINARY_DIR}/include/networkit)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
add_library(networkit
  external/networkit/networkit/cpp/algebraic/CSRMatrix.cpp
  external/networkit/networkit/cpp/algebraic/DenseMatrix.cpp
  external/networkit/networkit/cpp/algebraic/DynamicMatrix.cpp
  external/networkit/networkit/cpp/algebraic/Vector.cpp
  external/networkit/networkit/cpp/auxiliary/BloomFilter.cpp
  external/networkit/networkit/cpp/auxiliary/BucketPQ.cpp
  external/networkit/networkit/cpp/auxiliary/Log.cpp
  external/networkit/networkit/cpp/auxiliary/NumericTools.cpp
  external/networkit/networkit/cpp/auxiliary/Parallelism.cpp
  external/networkit/networkit/cpp/auxiliary/PrioQueueForInts.cpp
  external/networkit/networkit/cpp/auxiliary/Random.cpp
  external/networkit/networkit/cpp/auxiliary/SignalHandling.cpp
  external/networkit/networkit/cpp/auxiliary/Timer.cpp
  external/networkit/networkit/cpp/base/Algorithm.cpp
  external/networkit/networkit/cpp/centrality/ApproxBetweenness.cpp
  external/networkit/networkit/cpp/centrality/ApproxCloseness.cpp
  external/networkit/networkit/cpp/centrality/Betweenness.cpp
  external/networkit/networkit/cpp/centrality/Centrality.cpp
  external/networkit/networkit/cpp/centrality/Closeness.cpp
  external/networkit/networkit/cpp/centrality/CoreDecomposition.cpp
  external/networkit/networkit/cpp/centrality/DegreeCentrality.cpp
  external/networkit/networkit/cpp/centrality/DynApproxBetweenness.cpp
  external/networkit/networkit/cpp/centrality/DynBetweenness.cpp
  external/networkit/networkit/cpp/centrality/EigenvectorCentrality.cpp
  external/networkit/networkit/cpp/centrality/EstimateBetweenness.cpp
  external/networkit/networkit/cpp/centrality/KatzCentrality.cpp
  external/networkit/networkit/cpp/centrality/KPathCentrality.cpp
  external/networkit/networkit/cpp/centrality/LocalClusteringCoefficient.cpp
  external/networkit/networkit/cpp/centrality/LocalPartitionCoverage.cpp
  external/networkit/networkit/cpp/centrality/PageRank.cpp
  external/networkit/networkit/cpp/centrality/PermanenceCentrality.cpp
  external/networkit/networkit/cpp/centrality/Sfigality.cpp
  external/networkit/networkit/cpp/centrality/SpanningEdgeCentrality.cpp
  external/networkit/networkit/cpp/centrality/TopCloseness.cpp
  external/networkit/networkit/cpp/clique/MaxClique.cpp
  external/networkit/networkit/cpp/clique/MaximalCliques.cpp
  external/networkit/networkit/cpp/coarsening/ClusteringProjector.cpp
  external/networkit/networkit/cpp/coarsening/GraphCoarsening.cpp
  external/networkit/networkit/cpp/coarsening/MatchingCoarsening.cpp
  external/networkit/networkit/cpp/coarsening/ParallelPartitionCoarsening.cpp
  external/networkit/networkit/cpp/community/AdjustedRandMeasure.cpp
  external/networkit/networkit/cpp/community/ClusteringGenerator.cpp
  external/networkit/networkit/cpp/community/CommunityDetectionAlgorithm.cpp
  external/networkit/networkit/cpp/community/Conductance.cpp
  external/networkit/networkit/cpp/community/Coverage.cpp
  external/networkit/networkit/cpp/community/CoverHubDominance.cpp
  external/networkit/networkit/cpp/community/CutClustering.cpp
  external/networkit/networkit/cpp/community/DissimilarityMeasure.cpp
  external/networkit/networkit/cpp/community/DynamicNMIDistance.cpp
  external/networkit/networkit/cpp/community/EdgeCut.cpp
  external/networkit/networkit/cpp/community/GraphClusteringTools.cpp
  external/networkit/networkit/cpp/community/GraphStructuralRandMeasure.cpp
  external/networkit/networkit/cpp/community/HubDominance.cpp
  external/networkit/networkit/cpp/community/IntrapartitionDensity.cpp
  external/networkit/networkit/cpp/community/IsolatedInterpartitionConductance.cpp
  external/networkit/networkit/cpp/community/IsolatedInterpartitionExpansion.cpp
  external/networkit/networkit/cpp/community/JaccardMeasure.cpp
  external/networkit/networkit/cpp/community/LocalCoverEvaluation.cpp
  external/networkit/networkit/cpp/community/LocalPartitionEvaluation.cpp
  external/networkit/networkit/cpp/community/LPDegreeOrdered.cpp
  external/networkit/networkit/cpp/community/Modularity.cpp
  external/networkit/networkit/cpp/community/NMIDistance.cpp
  external/networkit/networkit/cpp/community/NodeStructuralRandMeasure.cpp
  external/networkit/networkit/cpp/community/ParallelAgglomerativeClusterer.cpp
  external/networkit/networkit/cpp/community/PartitionFragmentation.cpp
  external/networkit/networkit/cpp/community/PartitionHubDominance.cpp
  external/networkit/networkit/cpp/community/PartitionIntersection.cpp
  external/networkit/networkit/cpp/community/PLM.cpp
  external/networkit/networkit/cpp/community/PLP.cpp
  external/networkit/networkit/cpp/community/QualityMeasure.cpp
  external/networkit/networkit/cpp/community/SampledGraphStructuralRandMeasure.cpp
  external/networkit/networkit/cpp/community/SampledNodeStructuralRandMeasure.cpp
  external/networkit/networkit/cpp/community/StablePartitionNodes.cpp
  external/networkit/networkit/cpp/components/ConnectedComponents.cpp
  external/networkit/networkit/cpp/components/DynConnectedComponents.cpp
  external/networkit/networkit/cpp/components/ParallelConnectedComponents.cpp
  external/networkit/networkit/cpp/components/StronglyConnectedComponents.cpp
  external/networkit/networkit/cpp/components/WeaklyConnectedComponents.cpp
  external/networkit/networkit/cpp/correlation/Assortativity.cpp
  external/networkit/networkit/cpp/distance/AdamicAdarDistance.cpp
  external/networkit/networkit/cpp/distance/AlgebraicDistance.cpp
  external/networkit/networkit/cpp/distance/AllSimplePaths.cpp
  external/networkit/networkit/cpp/distance/APSP.cpp
  external/networkit/networkit/cpp/distance/BFS.cpp
  external/networkit/networkit/cpp/distance/CommuteTimeDistance.cpp
  external/networkit/networkit/cpp/distance/Diameter.cpp
  external/networkit/networkit/cpp/distance/Dijkstra.cpp
  external/networkit/networkit/cpp/distance/DynAPSP.cpp
  external/networkit/networkit/cpp/distance/DynBFS.cpp
  external/networkit/networkit/cpp/distance/DynDijkstra.cpp
  external/networkit/networkit/cpp/distance/DynSSSP.cpp
  external/networkit/networkit/cpp/distance/Eccentricity.cpp
  external/networkit/networkit/cpp/distance/EffectiveDiameter.cpp
  external/networkit/networkit/cpp/distance/EffectiveDiameterApproximation.cpp
  external/networkit/networkit/cpp/distance/GraphDistance.cpp
  external/networkit/networkit/cpp/distance/HopPlotApproximation.cpp
  external/networkit/networkit/cpp/distance/IncompleteDijkstra.cpp
  external/networkit/networkit/cpp/distance/JaccardDistance.cpp
  external/networkit/networkit/cpp/distance/NeighborhoodFunction.cpp
  external/networkit/networkit/cpp/distance/NeighborhoodFunctionApproximation.cpp
  external/networkit/networkit/cpp/distance/NeighborhoodFunctionHeuristic.cpp
  external/networkit/networkit/cpp/distance/NodeDistance.cpp
  external/networkit/networkit/cpp/distance/SSSP.cpp
  external/networkit/networkit/cpp/dynamics/DGSStreamParser.cpp
  external/networkit/networkit/cpp/dynamics/DGSWriter.cpp
  external/networkit/networkit/cpp/dynamics/GraphEvent.cpp
  external/networkit/networkit/cpp/dynamics/GraphEventHandler.cpp
  external/networkit/networkit/cpp/dynamics/GraphEventProxy.cpp
  external/networkit/networkit/cpp/dynamics/GraphUpdater.cpp
  external/networkit/networkit/cpp/edgescores/ChibaNishizekiQuadrangleEdgeScore.cpp
  external/networkit/networkit/cpp/edgescores/ChibaNishizekiTriangleEdgeScore.cpp
  external/networkit/networkit/cpp/edgescores/EdgeScore.cpp
  external/networkit/networkit/cpp/edgescores/EdgeScoreAsWeight.cpp
  external/networkit/networkit/cpp/edgescores/EdgeScoreBlender.cpp
  external/networkit/networkit/cpp/edgescores/EdgeScoreLinearizer.cpp
  external/networkit/networkit/cpp/edgescores/EdgeScoreNormalizer.cpp
  external/networkit/networkit/cpp/edgescores/GeometricMeanScore.cpp
  external/networkit/networkit/cpp/edgescores/PrefixJaccardScore.cpp
  external/networkit/networkit/cpp/edgescores/TriangleEdgeScore.cpp
  external/networkit/networkit/cpp/flow/EdmondsKarp.cpp
  external/networkit/networkit/cpp/generators/BarabasiAlbertGenerator.cpp
  external/networkit/networkit/cpp/generators/ChungLuGenerator.cpp
  external/networkit/networkit/cpp/generators/ClusteredRandomGraphGenerator.cpp
  external/networkit/networkit/cpp/generators/DorogovtsevMendesGenerator.cpp
  external/networkit/networkit/cpp/generators/DynamicBarabasiAlbertGenerator.cpp
  external/networkit/networkit/cpp/generators/DynamicDGSParser.cpp
  external/networkit/networkit/cpp/generators/DynamicDorogovtsevMendesGenerator.cpp
  external/networkit/networkit/cpp/generators/DynamicForestFireGenerator.cpp
  external/networkit/networkit/cpp/generators/DynamicGraphGenerator.cpp
  external/networkit/networkit/cpp/generators/DynamicGraphSource.cpp
  external/networkit/networkit/cpp/generators/DynamicHyperbolicGenerator.cpp
  external/networkit/networkit/cpp/generators/DynamicPathGenerator.cpp
  external/networkit/networkit/cpp/generators/DynamicPubWebGenerator.cpp
  external/networkit/networkit/cpp/generators/EdgeSwitchingMarkovChainGenerator.cpp
  external/networkit/networkit/cpp/generators/ErdosRenyiGenerator.cpp
  external/networkit/networkit/cpp/generators/HavelHakimiGenerator.cpp
  external/networkit/networkit/cpp/generators/HyperbolicGenerator.cpp
  external/networkit/networkit/cpp/generators/LFRGenerator.cpp
  external/networkit/networkit/cpp/generators/PowerlawDegreeSequence.cpp
  external/networkit/networkit/cpp/generators/PubWebGenerator.cpp
  external/networkit/networkit/cpp/generators/RegularRingLatticeGenerator.cpp
  external/networkit/networkit/cpp/generators/RmatGenerator.cpp
  external/networkit/networkit/cpp/generators/StaticDegreeSequenceGenerator.cpp
  external/networkit/networkit/cpp/generators/StaticGraphGenerator.cpp
  external/networkit/networkit/cpp/generators/StochasticBlockmodel.cpp
  external/networkit/networkit/cpp/generators/WattsStrogatzGenerator.cpp
  external/networkit/networkit/cpp/geometric/HyperbolicSpace.cpp
  external/networkit/networkit/cpp/global/ClusteringCoefficient.cpp
  external/networkit/networkit/cpp/global/GlobalClusteringCoefficient.cpp
  external/networkit/networkit/cpp/graph/Graph.cpp
  external/networkit/networkit/cpp/graph/GraphBuilder.cpp
  external/networkit/networkit/cpp/graph/GraphTools.cpp
  external/networkit/networkit/cpp/graph/KruskalMSF.cpp
  external/networkit/networkit/cpp/graph/RandomMaximumSpanningForest.cpp
  external/networkit/networkit/cpp/graph/RandomSpanningForest.cpp
  external/networkit/networkit/cpp/graph/Sampling.cpp
  external/networkit/networkit/cpp/graph/SpanningForest.cpp
  external/networkit/networkit/cpp/graph/UnionMaximumSpanningForest.cpp
  external/networkit/networkit/cpp/independentset/IndependentSetFinder.cpp
  external/networkit/networkit/cpp/independentset/Luby.cpp
  external/networkit/networkit/cpp/io/CoverReader.cpp
  external/networkit/networkit/cpp/io/CoverWriter.cpp
  external/networkit/networkit/cpp/io/DGSReader.cpp
  external/networkit/networkit/cpp/io/DibapGraphReader.cpp
  external/networkit/networkit/cpp/io/DotGraphWriter.cpp
  external/networkit/networkit/cpp/io/DotPartitionWriter.cpp
  external/networkit/networkit/cpp/io/DynamicGraphReader.cpp
  external/networkit/networkit/cpp/io/EdgeListCoverReader.cpp
  external/networkit/networkit/cpp/io/EdgeListPartitionReader.cpp
  external/networkit/networkit/cpp/io/EdgeListReader.cpp
  external/networkit/networkit/cpp/io/EdgeListWriter.cpp
  external/networkit/networkit/cpp/io/GMLGraphReader.cpp
  external/networkit/networkit/cpp/io/GMLGraphWriter.cpp
  external/networkit/networkit/cpp/io/GraphIO.cpp
  external/networkit/networkit/cpp/io/GraphToolBinaryReader.cpp
  external/networkit/networkit/cpp/io/GraphToolBinaryWriter.cpp
  external/networkit/networkit/cpp/io/KONECTGraphReader.cpp
  external/networkit/networkit/cpp/io/LineFileReader.cpp
  external/networkit/networkit/cpp/io/MatrixMarketReader.cpp
  external/networkit/networkit/cpp/io/METISGraphReader.cpp
  external/networkit/networkit/cpp/io/METISGraphWriter.cpp
  external/networkit/networkit/cpp/io/METISParser.cpp
  external/networkit/networkit/cpp/io/PartitionReader.cpp
  external/networkit/networkit/cpp/io/PartitionWriter.cpp
  external/networkit/networkit/cpp/io/RasterReader.cpp
  external/networkit/networkit/cpp/io/SNAPEdgeListPartitionReader.cpp
  external/networkit/networkit/cpp/io/SNAPGraphReader.cpp
  external/networkit/networkit/cpp/io/SNAPGraphWriter.cpp
  external/networkit/networkit/cpp/linkprediction/AdamicAdarIndex.cpp
  external/networkit/networkit/cpp/linkprediction/AdjustedRandIndex.cpp
  external/networkit/networkit/cpp/linkprediction/AlgebraicDistanceIndex.cpp
  external/networkit/networkit/cpp/linkprediction/CommonNeighborsIndex.cpp
  external/networkit/networkit/cpp/linkprediction/EvaluationMetric.cpp
  external/networkit/networkit/cpp/linkprediction/JaccardIndex.cpp
  external/networkit/networkit/cpp/linkprediction/KatzIndex.cpp
  external/networkit/networkit/cpp/linkprediction/LinkPredictor.cpp
  external/networkit/networkit/cpp/linkprediction/LinkThresholder.cpp
  external/networkit/networkit/cpp/linkprediction/MissingLinksFinder.cpp
  external/networkit/networkit/cpp/linkprediction/NeighborhoodDistanceIndex.cpp
  external/networkit/networkit/cpp/linkprediction/NeighborhoodUtility.cpp
  external/networkit/networkit/cpp/linkprediction/NeighborsMeasureIndex.cpp
  external/networkit/networkit/cpp/linkprediction/PrecisionRecallMetric.cpp
  external/networkit/networkit/cpp/linkprediction/PredictionsSorter.cpp
  external/networkit/networkit/cpp/linkprediction/PreferentialAttachmentIndex.cpp
  external/networkit/networkit/cpp/linkprediction/RandomLinkSampler.cpp
  external/networkit/networkit/cpp/linkprediction/ResourceAllocationIndex.cpp
  external/networkit/networkit/cpp/linkprediction/ROCMetric.cpp
  external/networkit/networkit/cpp/linkprediction/SameCommunityIndex.cpp
  external/networkit/networkit/cpp/linkprediction/TotalNeighborsIndex.cpp
  external/networkit/networkit/cpp/linkprediction/UDegreeIndex.cpp
  external/networkit/networkit/cpp/linkprediction/VDegreeIndex.cpp
  external/networkit/networkit/cpp/matching/LocalMaxMatcher.cpp
  external/networkit/networkit/cpp/matching/Matcher.cpp
  external/networkit/networkit/cpp/matching/Matching.cpp
  external/networkit/networkit/cpp/matching/PathGrowingMatcher.cpp
  external/networkit/networkit/cpp/numerics/LAMG/MultiLevelSetup.cpp
  external/networkit/networkit/cpp/numerics/LAMG/SolverLamg.cpp
  external/networkit/networkit/cpp/overlap/HashingOverlapper.cpp
  external/networkit/networkit/cpp/overlap/Overlapper.cpp
  external/networkit/networkit/cpp/scd/ApproximatePageRank.cpp
  external/networkit/networkit/cpp/scd/GCE.cpp
  external/networkit/networkit/cpp/scd/PageRankNibble.cpp
  external/networkit/networkit/cpp/scd/SelectiveCommunityDetector.cpp
  external/networkit/networkit/cpp/simulation/EpidemicSimulationSEIR.cpp
  external/networkit/networkit/cpp/sparsification/ChanceCorrectedTriangleScore.cpp
  external/networkit/networkit/cpp/sparsification/ForestFireScore.cpp
  external/networkit/networkit/cpp/sparsification/GlobalThresholdFilter.cpp
  external/networkit/networkit/cpp/sparsification/LocalDegreeScore.cpp
  external/networkit/networkit/cpp/sparsification/LocalSimilarityScore.cpp
  external/networkit/networkit/cpp/sparsification/MultiscaleScore.cpp
  external/networkit/networkit/cpp/sparsification/RandomEdgeScore.cpp
  external/networkit/networkit/cpp/sparsification/RandomNodeEdgeScore.cpp
  external/networkit/networkit/cpp/sparsification/SCANStructuralSimilarityScore.cpp
  external/networkit/networkit/cpp/sparsification/SimmelianOverlapScore.cpp
  external/networkit/networkit/cpp/sparsification/SimmelianScore.cpp
  external/networkit/networkit/cpp/sparsification/Sparsifiers.cpp
  external/networkit/networkit/cpp/structures/Cover.cpp
  external/networkit/networkit/cpp/structures/Partition.cpp
  external/networkit/networkit/cpp/structures/UnionFind.cpp
  external/networkit/networkit/cpp/viz/FruchtermanReingold.cpp
  external/networkit/networkit/cpp/viz/Layouter.cpp
  external/networkit/networkit/cpp/viz/MaxentStress.cpp
  external/networkit/networkit/cpp/viz/MultilevelLayouter.cpp
  external/networkit/networkit/cpp/viz/PivotMDS.cpp
  external/networkit/networkit/cpp/viz/PostscriptWriter.cpp
)

# Create library from common files to avoid re-compiling them 
add_library(seidr-common src/External/gzstream.C src/Common/common.cpp 
  src/Common/fs.cpp src/Common/mpims.cpp src/Common/BSlogger.cpp)
add_library(seidr-glmnet src/External/glmnet5.f90)

# Build executables:
## Seidr
add_executable(seidr src/Seidr/adjacency.cpp src/Seidr/seidr.cpp 
src/Seidr/viewRanks.cpp src/Seidr/import.cpp src/Seidr/aggregate.cpp
src/Seidr/threshold.cpp src/Seidr/roc.cpp  src/Seidr/convert.cpp 
src/Seidr/backbone.cpp src/Seidr/compare_clusters.cpp src/Seidr/Serialize.cpp
src/Seidr/stats.cpp src/Seidr/index.cpp src/Seidr/neighbours.cpp
src/Seidr/test.cpp src/Seidr/compare.cpp src/Seidr/resolve.cpp
src/Seidr/reheader.cpp)

## Tigress
add_executable(tigress src/Tigress/tiglm.cpp src/Tigress/tigress.cpp 
  src/Tigress/glmnet2.cpp)
## el-ensemble
add_executable(el-ensemble src/ElNet/el-ensemble.cpp src/ElNet/elnet-fun.cpp 
  src/ElNet/glmnetx.cpp)
## SVM ensemble
add_executable(svm-ensemble src/SVM/svm-ensemble.cpp src/SVM/svm-fun.cpp 
  src/External/svm.cpp)
## llr ensemble
add_executable(llr-ensemble src/SVM/linear-ensemble.cpp src/SVM/linear-fun.cpp 
  src/External/linear.cpp src/External/tron.cpp)
## Narromi
add_executable(narromi src/Narromi/stats_fun.cpp src/Narromi/IP_LPT.cpp 
  src/Narromi/narromi_fun.cpp src/Narromi/narromi.cpp)
## Correlation
add_executable(correlation src/Correlation/correlation.cpp 
  src/Correlation/cor_fun.cpp)
## Anoverence
add_executable(anoverence src/Anoverence/anoverence.cpp 
  src/Anoverence/anova-fun.cpp)
## CLR
add_executable(mi src/MI/mi.cpp src/MI/mi_fun.cpp)
## PLSNET
add_executable(plsnet src/PLSNET/plsnet.cpp src/PLSNET/plsnet-fun.cpp)
## Partial Correlation
add_executable(pcor src/PCor/pcor.cpp src/PCor/pcor-fun.cpp)
## Genie3
add_executable(genie3 src/GENIE3/genie3.cpp src/GENIE3/genie3-fun.cpp 
  external/ranger/src/Data.cpp external/ranger/src/DataDouble.cpp 
  external/ranger/src/DataChar.cpp external/ranger/src/DataFloat.cpp 
  external/ranger/src/Forest.cpp external/ranger/src/ForestRegression.cpp 
  external/ranger/src/Tree.cpp external/ranger/src/TreeRegression.cpp 
  external/ranger/src/utility.cpp)


# Link libraries
target_link_libraries(seidr seidr-common)
target_link_libraries(tigress seidr-common)
target_link_libraries(svm-ensemble seidr-common)
target_link_libraries(llr-ensemble seidr-common)
target_link_libraries(el-ensemble seidr-common)
target_link_libraries(narromi seidr-common)
target_link_libraries(correlation seidr-common)
target_link_libraries(anoverence seidr-common)
target_link_libraries(mi seidr-common)
target_link_libraries(plsnet seidr-common)
target_link_libraries(pcor seidr-common)
target_link_libraries(genie3 seidr-common)

target_link_libraries(tigress seidr-glmnet)
target_link_libraries(el-ensemble seidr-glmnet)

# Enable networkit and OpenMP
find_package(OpenMP COMPONENTS C CXX REQUIRED)
set_target_properties(seidr PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
set_target_properties(seidr PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
set_target_properties(networkit PROPERTIES COMPILE_FLAGS 
                      "${OpenMP_CXX_FLAGS} -DNDEBUG -DNOLOGGING")
set_target_properties(networkit PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
target_link_libraries(seidr networkit)
target_link_libraries(seidr ${OpenMP_CXX_LIBRARIES})

find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( seidr ${ZLIB_LIBRARIES} )

find_package( Armadillo REQUIRED )
include_directories( ${ARMADILLO_INCLUDE_DIRS})
if(NOT CMAKE_BUILD_STATIC)
  target_link_libraries( seidr ${ARMADILLO_LIBRARIES})
  target_link_libraries( tigress ${ARMADILLO_LIBRARIES})
  target_link_libraries( svm-ensemble ${ARMADILLO_LIBRARIES})
  target_link_libraries( llr-ensemble ${ARMADILLO_LIBRARIES})
  target_link_libraries( el-ensemble ${ARMADILLO_LIBRARIES})
  target_link_libraries( narromi ${ARMADILLO_LIBRARIES})
  target_link_libraries( correlation ${ARMADILLO_LIBRARIES})
  target_link_libraries( anoverence ${ARMADILLO_LIBRARIES})
  target_link_libraries( mi ${ARMADILLO_LIBRARIES})
  target_link_libraries( plsnet ${ARMADILLO_LIBRARIES})
  target_link_libraries( pcor ${ARMADILLO_LIBRARIES})
  target_link_libraries( genie3 ${ARMADILLO_LIBRARIES})
  find_package(BLAS REQUIRED)
  include_directories(${BLAS_INCLUDE_PATH})
  target_link_libraries(narromi ${BLAS_LIBRARIES})
  target_link_libraries(llr-ensemble ${BLAS_LIBRARIES})
else()
  find_package( BLAS REQUIRED )
  include_directories( ${BLAS_INCLUDE_DIRS})
  find_package( LAPACK REQUIRED )
  include_directories( ${LAPACK_INCLUDE_DIRS})
  target_link_libraries( seidr ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util hts z bz2 gsl lzma)
  target_link_libraries( tigress ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( svm-ensemble ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( llr-ensemble ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( el-ensemble ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( narromi ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( correlation ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( anoverence ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( mi ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( plsnet ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( pcor ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
  target_link_libraries( genie3 ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m quadmath mpi pthread dl rt open-pal open-rte util)
endif()

find_package( Boost COMPONENTS system filesystem REQUIRED )
include_directories( ${BOOST_INCLUDE_DIR} )
target_link_libraries( seidr ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries( tigress  ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries( svm-ensemble  ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries( llr-ensemble  ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries( el-ensemble  ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries( narromi ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries( correlation ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries( anoverence ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries( mi ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries( plsnet ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries( pcor ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries( genie3 ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
if(NOT CMAKE_BUILD_STATIC)
  target_link_libraries( seidr ${MPI_LIBRARIES})
  target_link_libraries( tigress ${MPI_LIBRARIES})
  target_link_libraries( svm-ensemble ${MPI_LIBRARIES})
  target_link_libraries( llr-ensemble ${MPI_LIBRARIES})
  target_link_libraries( el-ensemble ${MPI_LIBRARIES})
  target_link_libraries( narromi ${MPI_LIBRARIES})
  target_link_libraries( mi ${MPI_LIBRARIES})
  target_link_libraries( plsnet ${MPI_LIBRARIES})
  target_link_libraries( genie3 ${MPI_LIBRARIES})
  find_package(BLAS REQUIRED)
  include_directories(${BLAS_INCLUDE_PATH})
  target_link_libraries(narromi ${BLAS_LIBRARIES})
  target_link_libraries(llr-ensemble ${BLAS_LIBRARIES})

  find_package(HTSLIB REQUIRED)
  include_directories(${HTSLIB_INCLUDE_DIR})
  target_link_libraries(seidr ${HTSLIB_LIBRARIES})

  find_package(GSL REQUIRED)
  include_directories(${GSL_INCLUDE_PATH})
  target_link_libraries(seidr ${GSL_LIBRARIES})
endif()

if(MPI_COMPILE_FLAGS)
  set_target_properties(seidr PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(tigress PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(narromi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(svm-ensemble PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(llr-ensemble PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(el-ensemble PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(mi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(plsnet PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(genie3 PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(seidr PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(tigress PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(svm-ensemble PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(llr-ensemble PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(el-ensemble PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(narromi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(mi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(plsnet PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  set_target_properties(genie3 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(narromi glpk)
target_link_libraries(genie3 pthread)
