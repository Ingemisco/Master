add_compile_options(
  -O3 
	-march=native 
	-mtune=native
  -DNDEBUG
	-funroll-loops
  -fomit-frame-pointer
	# -fno-exceptions 
	# -fno-rtti
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-flto)
    add_link_options(-flto -fuse-ld=gold)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-flto=thin)
    add_link_options(-flto=thin -fuse-ld=lld)
endif()

set(DEBUG_VALUE 0)
file(GLOB SRC_FILES "src/*.cpp")
file(GLOB DATAGEN_SRC_FILES "datageneration/*.cpp" "src/datastructures.cpp" "src/distance.cpp")

add_executable(polyline ${SRC_FILES})
add_executable(datagen ${DATAGEN_SRC_FILES})

find_package(Boost REQUIRED COMPONENTS program_options)
target_link_libraries(polyline PRIVATE Boost::headers Boost::program_options)
target_link_libraries(datagen PRIVATE Boost::headers Boost::program_options)
target_link_libraries(polyline PUBLIC OpenMP::OpenMP_CXX)
