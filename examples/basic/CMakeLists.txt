cmake_minimum_required(VERSION 2.8)

project(basic-example C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

add_executable(basic-example main.c)
target_link_libraries(basic-example predict)
