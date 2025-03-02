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

set(DEBUG_VALUE 1)
file(GLOB SRC_FILES "src/*.cpp")
file(GLOB DATAGEN_SRC_FILES "datageneration/*.cpp" "src/datastructures.cpp" "src/distance.cpp")

add_executable(polyline ${SRC_FILES})
add_executable(datagen ${DATAGEN_SRC_FILES})

find_package(Boost REQUIRED COMPONENTS program_options)

target_link_libraries(polyline PRIVATE
  -fsanitize=address  
  -fsanitize=undefined
  Boost::headers 
  Boost::program_options
)

target_link_libraries(datagen PRIVATE
  -fsanitize=address  
  -fsanitize=undefined
  Boost::headers 
  Boost::program_options
)
