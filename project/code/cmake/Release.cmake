add_compile_options(
  -O3
  -DNDEBUG
)

file(GLOB SRC_FILES "src/*.cpp")
file(GLOB DATAGEN_SRC_FILES "datageneration/*.cpp")

add_executable(polyline ${SRC_FILES})
add_executable(datagen ${DATAGEN_SRC_FILES})

target_include_directories(polyline PRIVATE ${CMAKE_SOURCE_DIR}/include)

find_package(Boost REQUIRED COMPONENTS program_options)
target_link_libraries(polyline PRIVATE Boost::headers Boost::program_options)
target_link_libraries(datagen PRIVATE Boost::headers Boost::program_options)
