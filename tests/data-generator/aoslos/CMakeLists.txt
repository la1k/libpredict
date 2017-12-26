cmake_minimum_required(VERSION 2.8)

project(aoslos-data-generator CXX C)

set(LIBPREDICT_SOURCE_DIR "../../..")
include_directories("${LIBPREDICT_SOURCE_DIR}/tests" "${LIBPREDICT_SOURCE_DIR}/src")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

add_executable(aoslos-data-generator aoslos-data-generator.cpp ${LIBPREDICT_SOURCE_DIR}/test_utils/orbits_from_file.cpp ${LIBPREDICT_SOURCE_DIR}/src/pass_utils.c)
target_link_libraries(aoslos-data-generator predict)
