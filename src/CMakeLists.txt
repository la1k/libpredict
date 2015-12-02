add_library(predict SHARED orbit.c sun.c observer.c sdp4.c sgp4.c refraction.c unsorted.c julian_date.c version.c)

add_library(predict_static STATIC orbit.c sun.c observer.c sdp4.c sgp4.c refraction.c unsorted.c julian_date.c)

set_target_properties(predict_static PROPERTIES OUTPUT_NAME predict)

# The VERSION is the version number in the file name:
# libpredict.so.<VERSION>. This should be the full major.minor.patch
# version number.
set_target_properties(predict PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(predict_static PROPERTIES VERSION ${PROJECT_VERSION})

# The SOVERSION is the version number in the soname:
# libpredict.so.<SOVERSION>. This should be only the major version
# number, in order for applications to get a runtime dependency on an
# ABI version compatible with the library at build time.
set_target_properties(predict PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(predict_static PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})

# Header files needed to compile an application using libpredict
set(predict_public_header)
list(APPEND predict_public_header
  ${LIBPREDICT_INCLUDE_DIR}/predict/predict.h
)
set_target_properties(predict PROPERTIES PUBLIC_HEADER "${predict_public_header}")
# Maybe needed on Mac OS?
#set_target_properties(predict PROPERTIES FRAMEWORK TRUE)

target_link_libraries(predict m)
target_link_libraries(predict_static m)

# Tell the linker to use the symbol version script
get_property(predict_link_flags TARGET predict PROPERTY LINK_FLAGS)
set(predict_link_flags "${predict_link_flags} -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/libpredict.symver")
set_target_properties(predict PROPERTIES LINK_FLAGS ${predict_link_flags})


install(TARGETS predict predict_static
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/predict
)

# pkg-config file
configure_file(predict.pc.in predict.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/predict.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)


# python extension
option( PYTHON_BINDINGS "Build python bindings via SWIG" ON )
if( PYTHON_BINDINGS )
    find_package( SWIG )
    find_package( PythonLibs )
    find_package( PythonInterp )
endif ()

if( SWIG_FOUND AND PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND )
    include( UseSWIG )
    include_directories( ${PYTHON_INCLUDE_PATH} )
    include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
    swig_add_module( pypredict python pypredict.i julian_date.c observer.c orbit.c sdp4.c sgp4.c sun.c unsorted.c version.c )
    swig_link_libraries(pypredict ${PYTHON_LIBRARIES} )
    set_source_files_properties( pypredict.i PROPERTIES CPLUSPLUS ON )
    set_source_files_properties( pypredict.i PROPERTIES SWIG_FLAGS "-includeall" )

    execute_process( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print( sysconfig.get_python_lib( plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}' ) )"
        OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE )
    get_filename_component( _ABS_PYTHON_MODULE_PATH ${_ABS_PYTHON_MODULE_PATH} ABSOLUTE )
    file( RELATIVE_PATH _REL_PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX} ${_ABS_PYTHON_MODULE_PATH} )
    
    set( PYTHON_MODULE_PATH ${_REL_PYTHON_MODULE_PATH} )
    install( FILES ${CMAKE_CURRENT_BINARY_DIR}/_pypredict.so DESTINATION ${PYTHON_MODULE_PATH} )
    install( FILES ${CMAKE_CURRENT_BINARY_DIR}/pypredict.py DESTINATION ${PYTHON_MODULE_PATH} )
    
    set( LIBPREDICT_PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/${PYTHON_MODULE_PATH} )
    message( STATUS "Building python bindings via swig. Will be installed under ${LIBFTDI_PYTHON_MODULE_PATH}" )
else ()
    message(STATUS "Not building python bindings")
endif ()
