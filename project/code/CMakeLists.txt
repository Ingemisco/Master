cmake_minimum_required(VERSION 3.10)

project(polyline LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(INCLUDE_DIRECTORY include)

include_directories(${INCLUDE_DIRECTORY})
include_directories(${CMAKE_BINARY_DIR})

file(GLOB_RECURSE HEADERS ${INCLUDE_DIRECTORY}/*.h)
file(GLOB_RECURSE SOURCES src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Set compile flags for clangd
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_23)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30.0")
  cmake_policy(SET CMP0167 NEW)
endif()

set(Boost_USE_MULTITHREADED ON)
set(Boost_NO_WARN_NEW_VERSIONS 1)

find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)

target_link_libraries(${PROJECT_NAME} PRIVATE Boost::headers Boost::program_options )

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif()

target_compile_options(${PROJECT_NAME} PRIVATE
  -Wall
  -Wextra
  -Wpedantic
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wweak-vtables)
    if(WIN32)
        add_compile_options(-fms-extensions)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # using GCC
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(DEBUG_VALUE 1)

  target_compile_options(${PROJECT_NAME} PRIVATE
    -O0                       # No optimizations
    -g                        # Debug flags
    -fsanitize=address        # Enable AddressSanitizer
    -fsanitize=undefined      # Enable Undefined Behavior Sanitizer
    -fno-omit-frame-pointer   # Required for AddressSanitizer
  )

  target_link_libraries(${PROJECT_NAME} PRIVATE
    -fsanitize=address        # Link with AddressSanitizer
    -fsanitize=undefined      # Link with Undefined Behavior Sanitizer
  )

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wweak-vtables)
  endif ()

elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(DEBUG_VALUE 0)
  target_compile_options(${PROJECT_NAME} PRIVATE
    -O3 
    -DNDEBUG
  )
endif()

configure_file(
  ${INCLUDE_DIRECTORY}/config.h.in
  ${CMAKE_BINARY_DIR}/config.h
  @ONLY
)

