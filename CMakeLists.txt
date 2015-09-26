# The minimum required version hasn't been tested. Feel free to adjust
# downwards if necessary.
cmake_minimum_required(VERSION 2.8)

project(libpredict C CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x") #for tests written in C++

include(GNUInstallDirs)

# The project version number is also the version of the library ABI:
# major = bump on incompatible change
# minor = bump on compatible change
# patch = bump on any change

# These variables will be set automatically by newer CMake (at least
# 3.0) if the project command is given a VERSION argument. Set them
# manually for now.
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)
#PROJECT_VERSION_TWEAK isn't used
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
set(LIBPREDICT_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(LIBPREDICT_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(LIBPREDICT_VERSION_PATCH ${PROJECT_VERSION_MINOR})
#LIBPREDICT_VERSION_TWEAK isn't used
set(LIBPREDICT_VERSION ${PROJECT_VERSION})

set(LIBPREDICT_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(AFTER ${LIBPREDICT_INCLUDE_DIR})

# Tests and examples need to link with the library
set(LIBPREDICT_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/src)

# Testdata location
set(LIBPREDICT_TEST_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

add_subdirectory(include)
add_subdirectory(src)

# Add a target to generate API documentation with Doxygen
set(DOXYGEN_OUTPUT_DIR "doc")
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(Doxyfile.in Doxyfile @ONLY)
add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

add_subdirectory(travis)

enable_testing()
add_subdirectory(tests)
