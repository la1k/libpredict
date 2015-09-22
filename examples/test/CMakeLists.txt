cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_library(TestCase TestCase.cpp)

add_executable(TestSatellite TestSatellite.cpp)
target_link_libraries(TestSatellite TestCase predict)

add_executable(TestSun TestSun.cpp)
target_link_libraries(TestSun TestCase predict)

add_executable(TestMoon TestMoon.cpp)
target_link_libraries(TestMoon TestCase predict)
