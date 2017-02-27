cmake_minimum_required(VERSION 2.8)

project(passplot-example C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

add_executable(generate_pass_data generate_pass_data.c)
target_link_libraries(generate_pass_data predict m)
