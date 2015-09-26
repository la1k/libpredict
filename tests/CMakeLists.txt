link_directories(${LIBPREDICT_LIBRARY_DIR})

add_executable(dummy-t dummy-t.c)
add_test(NAME dummy COMMAND dummy-t)

add_executable(link-t link-t.c)
target_link_libraries(link-t predict)
add_test(NAME link COMMAND link-t)

#test against satellite data
add_executable(orbit-t orbit-t.cpp testcase_reader.cpp)
target_link_libraries(orbit-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/sat_*.test")
foreach(file ${files})
	get_filename_component(testname ${file} NAME_WE)
	add_test(NAME ${testname} COMMAND orbit-t ${file})
endforeach()

#test aos/los on satellite data
add_executable(aoslos-t aoslos-t.cpp testcase_reader.cpp)
target_link_libraries(aoslos-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/sat_*.test")
foreach(file ${files})
	get_filename_component(testname ${file} NAME_WE)
	add_test(NAME aoslos-${testname} COMMAND aoslos-t ${file})
endforeach()

#test against sun data
add_executable(sun-t sun-t.cpp testcase_reader.cpp)
target_link_libraries(sun-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/sun_*.test")
foreach(file ${files})
	get_filename_component(testname ${file} NAME_WE)
	add_test(NAME ${testname} COMMAND sun-t ${file})
endforeach()

#test against moon data
add_executable(moon-t moon-t.cpp testcase_reader.cpp)
target_link_libraries(moon-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/moon_*.test")
foreach(file ${files})
	get_filename_component(testname ${file} NAME_WE)
	add_test(NAME ${testname} COMMAND moon-t ${file})
endforeach()
