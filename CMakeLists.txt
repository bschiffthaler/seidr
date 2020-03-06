cmake_minimum_required(VERSION 3.9.4)
cmake_policy(SET CMP0074 NEW)

# Are we building a package? Link most things statically
if (CMAKE_LINK_STATIC)
  if (APPLE) # Use static libs from here
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/bundle/include")
    link_directories("${CMAKE_CURRENT_SOURCE_DIR}/bundle/lib")
    set(CMAKE_CXX_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/bundle/bin/clang++")
    set(CMAKE_C_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/bundle/bin/clang")
  endif()
  set(Boost_USE_STATIC_LIBS ON)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_USE_WRAPPER")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(HDF5_USE_STATIC_LIBRARIES ON)
  set(TBB_REQUIRE_STATIC ON)
  set(GLPK_REQUIRE_STATIC ON)
  set(CLP_REQUIRE_STATIC ON)
endif()

project(seidr)
include(ExternalProject)

enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

include_directories(src/Common src/Narromi src/Tigress src/SVM 
  src/ElNet src/Seidr src/Correlation src/MI src/PLSNET src/PCor 
  src/GENIE3 src/Anoverence external/bgzf)

#######################################
## Version, C++ Standard, build type ##
#######################################

# Define Version string
set(SEIDR_MAJOR_VERSION 0)
set(SEIDR_MINOR_VERSION 13)
set(SEIDR_PATCH_VERSION 0)
set(SEIDR_VERSION
  ${SEIDR_MAJOR_VERSION}.${SEIDR_MINOR_VERSION}.${SEIDR_PATCH_VERSION})

# C++ Standard
set(CMAKE_CXX_STANDARD 11)

# Explicitly set build type if user hasn't chosen one
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CPACK_PACKAGE_VENDOR "Bastian Schiffthaler")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Seidr")
set(CPACK_PACKAGE_VERSION_MAJOR ${SEIDR_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${SEIDR_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${SEIDR_PATCH_VERSION})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /dist
    /.*build.*
    /\\\\.DS_Store
)
include(CPack)

##########################
## Set up flags ##
##########################

set(CMAKE_OPTIMIZATION_LEVEL_PRIVATE "-O0")
if (CMAKE_BUILD_TYPE MATCHES "Release" OR CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
  set(CMAKE_OPTIMIZATION_LEVEL_PRIVATE "-O3")
endif()
if (CMAKE_OPTIMIZATION_LEVEL)
  set(CMAKE_OPTIMIZATION_LEVEL_PRIVATE "-O${CMAKE_OPTIMIZATION_LEVEL}")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE}")
set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE}")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE} -g")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-DDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE} -g")

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE} -g -Wall  -Wextra -Wpedantic")
set(CMAKE_Fortran_FLAGS_DEBUG "-DDEBUG ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE} -g -Wall  -Wextra -Wpedantic")

# Add common CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEIDR_VERSION=${SEIDR_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=${SEIDR_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_PRINT_ERRORS")

###################################################
## Compatibility with FORTRAN compiled libraries ##
###################################################

# Make sure FORTRAN has the same width in relevant datatypes to interface with 
# C++
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -integer-size 64 -autodouble -fixed")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8 -m64 -fdefault-real-8 -fdefault-double-8 -ffixed-form")
endif()

####################
## Compile faster ##
####################
if (CMAKE_PIPE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pipe")
endif()

##########################
## General Optimization ##
##########################

# Enable Processor tuning in command line (only in release mode)
if(CMAKE_BUILD_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -march=native -mtune=native")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -march=native -mtune=native")
endif()

# Allow user to select dangerous optimizations
if (CMAKE_OFAST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast")
endif()

# Can we use parallel sorting algorithms?
if (SEIDR_PSTL)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEIDR_PSTL")
endif()

###################
## Debug options ##
###################

# Enable extra debugging or inspection options
if(CMAKE_CLANG_ANALYZE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --analyze")
endif()

# Enable sanitizers
if(CMAKE_ASAN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
endif()

##########################
## Are we building MPI? ##
##########################
if(SEIDR_WITH_MPI)
  set(MPI_SRC "src/Common/mpiomp.cpp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEIDR_WITH_MPI")
else()
  set(MPI_SRC "src/Common/mpi_dummy.cpp")
endif()

# Enable building seidr with float as a default floating point datatype 
# in case memory is sparse
if(SEIDR_SCORE_FLOAT)
  ## Not supported anymore
  message(WARNING "SEIDR_SCORE_FLOAT is deprecated and was ignored")
endif()

# Enable OpenMP
find_package(OpenMP COMPONENTS CXX REQUIRED)
if(OPENMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Create library from common files to avoid re-compiling them 
add_library(seidr-common src/Common/common.cpp 
  src/Common/fs.cpp src/Common/BSlogger.cpp 
  "${MPI_SRC}")
add_library(seidr-glmnet external/glmnet/src/glmnet5dp.f)
add_library(seidr-bgzf external/bgzf/bgzf.c)

# Build a custom library from the networkit sources
add_subdirectory("${CMAKE_SOURCE_DIR}/deps")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/deps/include)

# Build executables:
## Seidr
add_executable(seidr src/Seidr/adjacency.cpp src/Seidr/aggregate.cpp 
               src/Seidr/asp.cpp src/Seidr/backbone.cpp 
               src/Seidr/compare.cpp src/Seidr/compare_clusters.cpp
               src/Seidr/convert.cpp src/Seidr/describe.cpp
               src/Seidr/graphstats.cpp src/Seidr/import.cpp
               src/Seidr/index.cpp src/Seidr/neighbours.cpp
               src/Common/parallel_control.cpp
               src/Seidr/reheader.cpp src/Seidr/resolve.cpp
               src/Seidr/roc.cpp src/Seidr/sample.cpp
               src/Seidr/Serialize.cpp src/Seidr/stats.cpp 
               src/Seidr/tau.cpp src/Seidr/test.cpp 
               src/Seidr/threshold.cpp src/Seidr/top.cpp
               src/Seidr/viewRanks.cpp src/Seidr/seidr.cpp)
## Tigress
add_executable(tigress src/Tigress/tiglm.cpp src/Tigress/tigress.cpp 
  src/Tigress/glmnet2.cpp)
## el-ensemble
add_executable(el-ensemble src/ElNet/el-ensemble.cpp src/ElNet/elnet-fun.cpp 
  src/ElNet/glmnetx.cpp)
## SVM ensemble
add_executable(svm-ensemble src/SVM/svm-ensemble.cpp src/SVM/svm-fun.cpp 
  external/libsvm/svm.cpp)
target_include_directories(svm-ensemble PRIVATE external/libsvm)
## llr ensemble
add_executable(llr-ensemble src/SVM/linear-ensemble.cpp src/SVM/linear-fun.cpp 
  external/liblinear/linear.cpp external/liblinear/tron.cpp)
target_include_directories(llr-ensemble PRIVATE external/liblinear)
## Narromi
if (NARROMI_USE_CLP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNARROMI_USE_CLP")
  set(IP_SRC "src/Narromi/IP_LPT2.cpp")
else()
  set(IP_SRC "src/Narromi/IP_LPT.cpp")
endif()
add_executable(narromi src/Narromi/stats_fun.cpp ${IP_SRC}
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
target_include_directories(genie3 PRIVATE external/ranger/src)

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

# Packages:

# ZLIB & BGZIP compression
target_include_directories( seidr PRIVATE external/bgzf )
target_link_libraries( seidr seidr-bgzf )
# Link zlib _after_ bgzf to avoid undef refl issues
find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( seidr ${ZLIB_LIBRARIES} )

if (SEIDR_PSTL)
  find_package(TBB REQUIRED)
  target_include_directories(seidr PRIVATE ${TBB_INCLUDE_DIR})
  target_link_libraries(seidr ${TBB_LIBRARIES})
endif()

if (CMAKE_LINK_STATIC)
   # Armadillo
  find_package( BLAS REQUIRED )
  include_directories( ${BLAS_INCLUDE_DIRS})
  target_link_libraries( seidr ${BLAS_LIBRARIES})
  target_link_libraries( tigress ${BLAS_LIBRARIES})
  target_link_libraries( svm-ensemble ${BLAS_LIBRARIES})
  target_link_libraries( llr-ensemble ${BLAS_LIBRARIES})
  target_link_libraries( el-ensemble ${BLAS_LIBRARIES})
  target_link_libraries( narromi ${BLAS_LIBRARIES})
  target_link_libraries( correlation ${BLAS_LIBRARIES})
  target_link_libraries( anoverence ${BLAS_LIBRARIES})
  target_link_libraries( mi ${BLAS_LIBRARIES})
  target_link_libraries( plsnet ${BLAS_LIBRARIES})
  target_link_libraries( pcor ${BLAS_LIBRARIES})
  target_link_libraries( genie3 ${BLAS_LIBRARIES})
else()
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
endif()

if (CMAKE_LINK_STATIC)
  find_package(HDF5 REQUIRED)
  include_directories( ${HDF5_INCLUDE_DIRS})
  target_link_libraries( seidr ${HDF5_LIBRARIES})
  target_link_libraries( tigress ${HDF5_LIBRARIES})
  target_link_libraries( svm-ensemble ${HDF5_LIBRARIES})
  target_link_libraries( llr-ensemble ${HDF5_LIBRARIES})
  target_link_libraries( el-ensemble ${HDF5_LIBRARIES})
  target_link_libraries( narromi ${HDF5_LIBRARIES})
  target_link_libraries( correlation ${HDF5_LIBRARIES})
  target_link_libraries( anoverence ${HDF5_LIBRARIES})
  target_link_libraries( mi ${HDF5_LIBRARIES})
  target_link_libraries( plsnet ${HDF5_LIBRARIES})
  target_link_libraries( pcor ${HDF5_LIBRARIES})
  target_link_libraries( genie3 ${HDF5_LIBRARIES})
endif()

# Boost
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_RELEASE_LIBS ON)
find_package( Boost COMPONENTS system filesystem program_options serialization REQUIRED )
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

if(SEIDR_WITH_MPI)
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
endif()


# Need to link BLAS because of liblinear/libsvm
find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_PATH})
target_link_libraries(llr-ensemble ${BLAS_LIBRARIES})
target_link_libraries(svm-ensemble ${BLAS_LIBRARIES})


if (NARROMI_USE_CLP)
  find_package(CLP REQUIRED)
  include_directories(${CLP_INCLUDE_PATH})
  target_link_libraries(narromi ${CLP_LIBRARIES})
  if (CMAKE_LINK_STATIC)
    find_package(BZip2 REQUIRED)
    target_link_libraries(narromi ${BZIP2_LIBRARIES})
    message(STATUS "BZ2: ${BZIP2_LIBRARIES}")
  endif()
else()
  # GLPK
  find_package(GLPK REQUIRED)
  include_directories(${GLPK_INCLUDE_PATH})
  target_link_libraries(narromi ${GLPK_LIBRARIES})
endif()

#####
## Install
#####

install(TARGETS anoverence correlation el-ensemble genie3 llr-ensemble mi
        narromi pcor plsnet seidr svm-ensemble tigress
        DESTINATION "bin")

message(STATUS "Summary:")
message(STATUS "  Seidr version: ${SEIDR_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
if (SEIDR_WITH_MPI)
  message(STATUS "  MPI: Yes")
else()
  message(STATUS "  MPI: No")
endif()
if (SEIDR_PSTL)
  message(STATUS "  PSTL: Yes")
else()
  message(STATUS "  PSTL: No")
endif()
if (CMAKE_ASAN)
  message(STATUS "  ASAN: Yes")
else()
  message(STATUS "  ASAN: No")
endif()
if (CMAKE_CLANG_ANALYZE)
  message(STATUS "  Clang Analyze: Yes")
else()
  message(STATUS "  Clang Analyze: No")
endif()
message(STATUS "  Optimization: ${CMAKE_OPTIMIZATION_LEVEL_PRIVATE}")
if (CMAKE_BUILD_NATIVE)
  message(STATUS "  Native build: Yes")
else()
  message(STATUS "  Native build: No")
endif()
if (NARROMI_USE_CLP)
  message(STATUS "  LP Engine: CLP")
  message(STATUS "    Linking: ${CLP_LIBRARIES}")
else()
  message(STATUS "  LP Engine: GLPK")
  message(STATUS "    Linking: ${GLPK_LIBRARIES}")
endif()
message(STATUS "  Armadillo: ${ARMADILLO_LIBRARIES}")
message(STATUS "  BLAS: ${BLAS_LIBRARIES}")
message(STATUS "  Boost: ${Boost_LIBRARIES}")
message(STATUS "  ZLIB: ${ZLIB_LIBRARIES}")
if (SEIDR_PSTL)
  message(STATUS "  TBB: ${TBB_LIBRARIES}")
endif()
if (CMAKE_BUILD_TYPE MATCHES "Release")
  message(STATUS "  CXX Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  message(STATUS "  Fortran flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
  message(STATUS "  CXX Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  message(STATUS "  Fortran flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
else()
  message(STATUS "  CXX Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  message(STATUS "  Fortran flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
endif()
