# The minimum required version hasn't been tested. Feel free to adjust
# downwards if necessary.
cmake_minimum_required(VERSION 2.8)

project(libpredict C)

# The project version number is also the version of the library ABI:
# major = bump on incompatible change
# minor = bump on compatible change
# patch = bump on any change

# These variables will be set automatically by newer CMake (at least
# 3.0) if the project command is given a VERSION argument. Set them
# manually for now.
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_PATCH 0)
#PROJECT_VERSION_TWEAK isn't used
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
set(LIBPREDICT_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(LIBPREDICT_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(LIBPREDICT_VERSION_PATCH ${PROJECT_VERSION_MINOR})
#LIBPREDICT_VERSION_TWEAK isn't used
set(LIBPREDICT_VERSION ${PROJECT_VERSION})

set(LIBPREDICT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(AFTER ${LIBPREDICT_INCLUDE_DIR})

add_subdirectory(src)
