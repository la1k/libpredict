cmake_minimum_required(VERSION 2.8)

project(setrise-example C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

add_executable(setrise-example main.c)
target_link_libraries(setrise-example predict)
