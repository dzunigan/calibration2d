project(calibration2d CXX)

cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

## Set CXX optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fext-numeric-literals -Wpedantic -Wall -Wextra")
## -mtune=native -march=native

## Set build type to Release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif(CMAKE_BUILD_TYPE STREQUAL "Release")
message(STATUS "Flags: " ${CMAKE_CXX_FLAGS})

## Build shared libs by default
set(BUILD_SHARED_LIBS "ON" CACHE BOOL "Build libraries as shared?")

## Dependencies options
set(GFLAGS_PREFER_EXPORTED_GFLAGS_CMAKE_CONFIGURATION OFF CACHE BOOL "Use exported cmake config to find gflags?")

## Dependencies
find_package(Boost REQUIRED filesystem system)
find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Gflags REQUIRED)
find_package(Glog REQUIRED)

option(BUILD_SYNC_TOOL "Build sync tool" ON)
if(BUILD_SYNC_TOOL)
  add_subdirectory(tools/sync)
endif(BUILD_SYNC_TOOL)

include_directories(include ${Boost_INCLUDE_DIRS} ${CERES_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GFLAGS_INCLUDE_DIR} ${GLOG_INCLUDE_DIRS})

## Header files
file(GLOB_RECURSE HEADER_FILES include/*.hpp include/*.h)
add_custom_target(calibrate_header_files SOURCES ${HEADER_FILES})

add_executable(calibrate src/calibrate.cpp src/random.cpp src/random_sampler.cpp src/support_measurement.cpp src/estimators.cpp)
target_link_libraries(calibrate ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${GLOG_LIBRARIES})

