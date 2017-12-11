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

#test against venus data
add_executable(venus-t venus-t.cpp testcase_reader.cpp)
target_link_libraries(venus-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/venus_*.test")
foreach(file ${files})
	get_filename_component(testname ${file} NAME_WE)
	add_test(NAME ${testname} COMMAND venus-t ${file})
endforeach()

#test against jupiter data
add_executable(jupiter-t jupiter-t.cpp testcase_reader.cpp)
target_link_libraries(jupiter-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/jupiter_*.test")
foreach(file ${files})
	get_filename_component(testname ${file} NAME_WE)
	add_test(NAME ${testname} COMMAND jupiter-t ${file})
endforeach()

#test max elevation function
add_executable(maxelevation-t maxelevation-t.cpp testcase_reader.cpp)
target_link_libraries(maxelevation-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/sat_*.test")
foreach(file ${files})
	get_filename_component(testname ${file} NAME_WE)
	add_test(NAME maxelevation-${testname} COMMAND maxelevation-t ${file})
endforeach()

#test predict_is_geostationary
add_executable(geostationary-t geostationary-t.cpp testcase_reader.cpp)
target_link_libraries(geostationary-t predict)
file(GLOB files "${LIBPREDICT_TEST_DATA_DIR}/sat_*.test")
add_test(NAME geostationary-celestrak-geostationary COMMAND geostationary-t ${LIBPREDICT_TEST_DATA_DIR}/geostationary.tle)
add_test(NAME geostationary-all-tles COMMAND geostationary-t ${LIBPREDICT_TEST_DATA_DIR}/large-tle-collection.tle ${LIBPREDICT_TEST_DATA_DIR}/geosynchronous_satellite_numbers.dat)
foreach(file ${files})
	get_filename_component(testname ${file} NAME_WE)
endforeach()
