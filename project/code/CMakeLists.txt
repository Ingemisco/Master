cmake_minimum_required(VERSION 3.14)
project(PolylineSimplification LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(DEBUG_VALUE 1)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30.0")
  cmake_policy(SET CMP0167 NEW)
endif()

# common compile flags
add_compile_options(
  -Wall 
  -Wextra 
  -Wpedantic
)

# include directories 
include_directories(${CMAKE_SOURCE_DIR})

# config file 
configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_SOURCE_DIR}/config.h @ONLY)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  include(${CMAKE_SOURCE_DIR}/cmake/Debug.cmake)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  include(${CMAKE_SOURCE_DIR}/cmake/Release.cmake)
elseif(CMAKE_BUILD_TYPE STREQUAL "Test")
  include(${CMAKE_SOURCE_DIR}/cmake/Test.cmake)
else()
  message(FATAL_ERROR "Unknown build type: ${CMAKE_BUILD_TYPE}")
endif()

