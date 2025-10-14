set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30.0")
  cmake_policy(SET CMP0167 NEW)
endif()

file(GLOB TEST_SRC_FILES *.cpp)
file(GLOB_RECURSE SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*/main\\.cpp$")

file(GLOB_RECURSE DATA_SOURCE_FILES ${CMAKE_SOURCE_DIR}/datageneration/*.cpp)
list(FILTER DATA_SOURCE_FILES EXCLUDE REGEX ".*/main\\.cpp$")

file(COPY ${CMAKE_SOURCE_DIR}/test_data/ DESTINATION ${CMAKE_BINARY_DIR}/tests/test_data)


add_library(polyline_lib STATIC ${SOURCE_FILES})
add_library(data_gen STATIC ${DATA_SOURCE_FILES})
target_include_directories(polyline_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(data_gen PUBLIC ${CMAKE_SOURCE_DIR}/datageneration)
target_link_libraries(polyline_lib PUBLIC OpenMP::OpenMP_CXX)

foreach(test_file ${TEST_SRC_FILES})
	get_filename_component(test_name ${test_file} NAME_WE)

	add_executable(${test_name} ${test_file})

	target_include_directories(${test_name} PRIVATE
		$(CMAKE_SOURCE_DIR)/include 
	)

	target_link_libraries(${test_name} PRIVATE 
		-fsanitize=address  
		-fsanitize=undefined
		Boost::headers 
		Boost::unit_test_framework
		polyline_lib
		data_gen
	)

	add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
