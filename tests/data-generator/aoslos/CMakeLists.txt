cmake_minimum_required(VERSION 2.8)

project(aoslos-data-generator CXX)

include_directories("../../" "../../../")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

add_executable(aoslos-data-generator aoslos-data-generator.cpp ../../test_utils/orbits_from_file.cpp)
target_link_libraries(aoslos-data-generator predict)
