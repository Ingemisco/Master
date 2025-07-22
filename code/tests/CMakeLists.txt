file(GLOB TEST_SRC_FILES *.cpp)
file(GLOB_RECURSE SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)

file(COPY ${CMAKE_SOURCE_DIR}/test_data/ DESTINATION ${CMAKE_BINARY_DIR}/tests/test_data)
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*/main\\.cpp$")

add_library(polyline_lib STATIC ${SOURCE_FILES})
target_include_directories(polyline_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)
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
	)

	add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
