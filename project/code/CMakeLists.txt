cmake_minimum_required(VERSION 3.10)

project(polyline LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories(include)

file(GLOB_RECURSE HEADERS include/*.h)

set(MAIN_CPP_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

file(GLOB_RECURSE SOURCES src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Set compile flags for clangd
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile flags for clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.30.0")
    cmake_policy(SET CMP0167 NEW)
endif ()

set(Boost_USE_MULTITHREADED ON)
set(Boost_NO_WARN_NEW_VERSIONS 1)

find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)

target_link_libraries(${PROJECT_NAME} PRIVATE Boost::headers Boost::program_options )

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${PROJECT_NAME} PRIVATE
        -Wall                     # Enable all warnings
        -Wextra                   # Enable extra warnings
        -Wpedantic                # Enable pedantic warnings
        -fsanitize=address        # Enable AddressSanitizer
        -fsanitize=undefined      # Enable Undefined Behavior Sanitizer
        -fno-omit-frame-pointer   # Required for AddressSanitizer
        -frtti
    )

    # Link the sanitizers
    target_link_libraries(${PROJECT_NAME} PRIVATE
        -fsanitize=address        # Link with AddressSanitizer
        -fsanitize=undefined      # Link with Undefined Behavior Sanitizer
    )
endif()
