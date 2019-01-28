cmake_minimum_required(VERSION 3.9.4)
project(seidr)
include(ExternalProject)

enable_language(Fortran)
link_directories(/mnt/picea/home/bastian/lib)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

include_directories(src/External src/Common src/Narromi src/Tigress src/SVM 
  src/ElNet src/Seidr src/Correlation src/MI src/PLSNET src/PCor 
  external/ranger/src src/GENIE3 src/Anoverence)

# Define Version string
set(SEIDR_MAJOR_VERSION 0)
set(SEIDR_MINOR_VERSION 10)
set(SEIDR_PATCH_VERSION 0)
set(SEIDR_VERSION
  ${SEIDR_MAJOR_VERSION}.${SEIDR_MINOR_VERSION}.${SEIDR_PATCH_VERSION})

# Add common CXX flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEIDR_VERSION=${SEIDR_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=${SEIDR_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_32BIT_WORD -DARMA_DONT_PRINT_ERRORS")

# Standard
set(CMAKE_CXX_STANDARD 11)

# Can we use parallel sorting algorithms?
if (CMAKE_PARALLEL_SORT) # No free compiler support except intel and msvc
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEIDR_PARALLEL_SORT")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# Default flags
if (CMAKE_OFAST)
  set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

# Make sure FORTRAN has the same width in relevant datatypes to interface with 
# C++
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-m64 -autodouble -fixed")
else()
  set(CMAKE_Fortran_FLAGS "-m64 -fdefault-real-8 -fdefault-double-8 -ffixed-form")
endif()

# Make sure we can use ASAN and undefined sanitizers in debug builds
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
if(CMAKE_CLANG_ANALYZE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --analyze")
endif()

# Enable Processor tuning in command line
if(CMAKE_BUILD_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  set(CMAKE_FORTRAN_FLAGS_RELEASE "${CMAKE_FORTRAN_FLAGS_RELEASE} -march=native")
endif()

# Enable building seidr with float as a default floating point datatype 
# in case memory is sparse
if(NOT SEIDR_SCORE_FLOAT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEIDR_SCORE_DOUBLE")
endif()

# Enable OpenMP
find_package(OpenMP COMPONENTS C CXX REQUIRED)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Create library from common files to avoid re-compiling them 
add_library(seidr-common src/External/gzstream.C src/Common/common.cpp 
  src/Common/fs.cpp src/Common/mpims.cpp src/Common/BSlogger.cpp)
add_library(seidr-glmnet src/External/glmnet5.f90)

# Build a custom library from the networkit sources
add_subdirectory("${CMAKE_SOURCE_DIR}/deps")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/deps/include)

# Build executables:
## Seidr
add_executable(seidr src/Seidr/adjacency.cpp src/Seidr/backbone.cpp 
src/Seidr/seidr.cpp src/Seidr/graphstats.cpp 
src/Seidr/viewRanks.cpp src/Seidr/import.cpp src/Seidr/aggregate.cpp
src/Seidr/threshold.cpp src/Seidr/roc.cpp  src/Seidr/convert.cpp 
src/Seidr/compare_clusters.cpp src/Seidr/Serialize.cpp
src/Seidr/stats.cpp src/Seidr/index.cpp src/Seidr/neighbours.cpp
src/Seidr/test.cpp src/Seidr/compare.cpp src/Seidr/resolve.cpp
src/Seidr/reheader.cpp src/Seidr/top.cpp src/Seidr/tau.cpp src/Seidr/asp.cpp
src/Seidr/sample.cpp)
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
add_executable(genie3 src/GENIE3/genie3.cpp src/GENIE3/genie3-fun.cpp)

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

target_link_libraries(seidr networkit)

target_link_libraries(seidr bs)

target_link_libraries(genie3 ranger pthread)

# Make libraries position independent
set_target_properties(seidr-common PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -fPIC")
set_target_properties(seidr-glmnet PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -fPIC") 

# Link TBB for intel parallel
if (CMAKE_PARALLEL_SORT)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_link_libraries(seidr tbb)
  endif()
endif()                                         

# Packages:

# ZLIB for BGZIP compression
find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( seidr ${ZLIB_LIBRARIES} )

# Armadillo
find_package( Armadillo REQUIRED )
include_directories( ${ARMADILLO_INCLUDE_DIRS})
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

# Boost
find_package( Boost COMPONENTS system filesystem program_options REQUIRED )
include_directories( ${BOOST_INCLUDE_DIR} )
target_link_libraries( seidr ${Boost_LIBRARIES})
target_link_libraries( tigress ${Boost_LIBRARIES})
target_link_libraries( svm-ensemble ${Boost_LIBRARIES})
target_link_libraries( llr-ensemble ${Boost_LIBRARIES})
target_link_libraries( el-ensemble ${Boost_LIBRARIES})
target_link_libraries( narromi ${Boost_LIBRARIES})
target_link_libraries( correlation ${Boost_LIBRARIES})
target_link_libraries( anoverence ${Boost_LIBRARIES})
target_link_libraries( mi ${Boost_LIBRARIES})
target_link_libraries( plsnet ${Boost_LIBRARIES})
target_link_libraries( pcor ${Boost_LIBRARIES})
target_link_libraries( genie3 ${Boost_LIBRARIES})

# OpenMPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries( seidr ${MPI_LIBRARIES})
target_link_libraries( tigress ${MPI_LIBRARIES})
target_link_libraries( svm-ensemble ${MPI_LIBRARIES})
target_link_libraries( llr-ensemble ${MPI_LIBRARIES})
target_link_libraries( el-ensemble ${MPI_LIBRARIES})
target_link_libraries( narromi ${MPI_LIBRARIES})
target_link_libraries( mi ${MPI_LIBRARIES})
target_link_libraries( plsnet ${MPI_LIBRARIES})
target_link_libraries( genie3 ${MPI_LIBRARIES})
target_link_libraries( anoverence ${MPI_LIBRARIES})
target_link_libraries( correlation ${MPI_LIBRARIES})
target_link_libraries( pcor ${MPI_LIBRARIES})

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

# BLAS (check if required?)
find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_PATH})
target_link_libraries(narromi ${BLAS_LIBRARIES})
target_link_libraries(llr-ensemble ${BLAS_LIBRARIES})

# HTSLIB for BGZIP
find_package(HTSLIB REQUIRED)
include_directories(${HTSLIB_INCLUDE_DIR})
target_link_libraries(seidr ${HTSLIB_LIBRARIES})

# GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_PATH})
target_link_libraries(seidr ${GSL_LIBRARIES})

find_package(GLPK REQUIRED)
include_directories(${GLPK_INCLUDE_PATH})
target_link_libraries(narromi ${GLPK_LIBRARIES})

target_link_libraries( seidr tcmalloc)
target_link_libraries( tigress tcmalloc)
target_link_libraries( svm-ensemble tcmalloc)
target_link_libraries( llr-ensemble tcmalloc)
target_link_libraries( el-ensemble tcmalloc)
target_link_libraries( narromi tcmalloc)
target_link_libraries( mi tcmalloc)
target_link_libraries( plsnet tcmalloc)
target_link_libraries( genie3 tcmalloc)
target_link_libraries( anoverence tcmalloc)
target_link_libraries( correlation tcmalloc)
target_link_libraries( pcor tcmalloc)
