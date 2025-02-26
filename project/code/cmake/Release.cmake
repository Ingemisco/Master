add_compile_options(
  -O3
  -DNDEBUG
)

set(DEBUG_VALUE 0)
file(GLOB SRC_FILES "src/*.cpp")
file(GLOB DATAGEN_SRC_FILES "datageneration/*.cpp")

add_executable(polyline ${SRC_FILES})
add_executable(datagen ${DATAGEN_SRC_FILES})

find_package(Boost REQUIRED COMPONENTS program_options)
target_link_libraries(polyline PRIVATE Boost::headers Boost::program_options)
target_link_libraries(datagen PRIVATE Boost::headers Boost::program_options)
