# The minimum required version hasn't been tested. Feel free to adjust
# downwards if necessary.
cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0.0)
# Set VERSION variables based on project().
cmake_policy(SET CMP0048 NEW)

# The project version number is also the version of the library ABI:
# major.minor.patch(.tweak)
# major = bump on incompatible change
# minor = bump on compatible change
# patch = bump on any change
# (tweak is not used by libpredict)
project(libpredict VERSION 3.0.0 LANGUAGES C CXX)
add_definitions(-D_XOPEN_SOURCE=700)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x") #for tests written in C++

include(GNUInstallDirs)
include(CheckCXXSourceCompiles)

set(LIBPREDICT_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(AFTER ${LIBPREDICT_INCLUDE_DIR})

# Tests and examples need to link with the library
set(LIBPREDICT_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/src)

# Testdata location
set(LIBPREDICT_TEST_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

set(CMAKE_REQUIRED_FLAGS "-Wno-unused-command-line-argument -Wl,--version-script ${CMAKE_CURRENT_SOURCE_DIR}/has_linker_version_script.symver")
check_cxx_source_compiles("int main(){return 0;}" HAS_LINKER_VERSION_SCRIPT)
set(CMAKE_REQUIRED_FLAGS "-Wno-unused-command-line-argument -Wl,-exported_symbols_list ${CMAKE_CURRENT_SOURCE_DIR}/has_linker_exported_symbols_list.symbols")
check_cxx_source_compiles("int main(){return 0;}" HAS_LINKER_EXPORTED_SYMBOLS_LIST)

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
