add_compile_options(
  -O0 
  -g 
  -fsanitize=address 
  -fsanitize=undefined 
  -fno-omit-frame-pointer
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wweak-vtables)
endif()

file(GLOB TEST_SRC_FILES "test/*.cpp")

target_include_directories(polyline PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_executable(test_algorithms ${TEST_SRC_FILES})

find_package(Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries(test_algorithms PRIVATE 
  -fsanitize=address  
  -fsanitize=undefined
  Boost::headers 
  Boost::unit_test_framework
)

