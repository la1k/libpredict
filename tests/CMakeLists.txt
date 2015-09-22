link_directories(${LIBPREDICT_LIBRARY_DIR})

add_executable(dummy-t dummy-t.c)
add_test(NAME dummy COMMAND dummy-t)

add_executable(link-t link-t.c)
target_link_libraries(link-t predict m)
add_test(NAME link COMMAND link-t)

#test against satellite data
add_executable(orbit-t orbit-t.cpp testcase.cpp)
target_compile_options(orbit-t PRIVATE "-std=c++11")
target_link_libraries(orbit-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/testcase_satellite_*")
foreach(file ${files})
	add_test(NAME orbit-${file} COMMAND orbit-t ${file})
endforeach()

#test against sun data
add_executable(sun-t sun-t.cpp testcase.cpp)
target_compile_options(sun-t PRIVATE "-std=c++11")
target_link_libraries(sun-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/testcase_sun_*")
foreach(file ${files})
	add_test(NAME sun-${file} COMMAND sun-t ${file})
endforeach()

#test against moon data
add_executable(moon-t moon-t.cpp testcase.cpp)
target_compile_options(moon-t PRIVATE "-std=c++11")
target_link_libraries(moon-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/testcase_moon_*")
foreach(file ${files})
	add_test(NAME moon-${file} COMMAND moon-t ${file})
endforeach()
