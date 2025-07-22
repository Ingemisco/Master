message(STATUS "Test mode enabled")

set(DEBUG_VALUE 0)
add_compile_options(
  -O3 
  -DNDEBUG
)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

enable_testing()

add_subdirectory(${CMAKE_SOURCE_DIR}/tests)


