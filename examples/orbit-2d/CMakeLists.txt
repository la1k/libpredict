cmake_minimum_required(VERSION 2.8)

project(orbit-plot-example CXX)

add_executable(generate_orbit_data generate_orbit_data.cpp)
target_link_libraries(generate_orbit_data predict)
