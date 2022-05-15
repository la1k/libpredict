add_library(predict SHARED orbit.c sun.c observer.c sdp4.c sgp4.c refraction.c unsorted.c julian_date.c version.c moon.c)

add_library(predict_static STATIC orbit.c sun.c observer.c sdp4.c sgp4.c refraction.c unsorted.c julian_date.c moon.c)

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

# Tell the linker to use the symbol version script, or export list if symbol versions aren't supported
get_property(predict_link_flags TARGET predict PROPERTY LINK_FLAGS)
if (${HAS_LINKER_VERSION_SCRIPT})
  set(predict_link_flags "${predict_link_flags} -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/libpredict.symver")
elseif(${HAS_LINKER_EXPORTED_SYMBOLS_LIST})
  set(predict_link_flags "${predict_link_flags} -Wl,-exported_symbols_list ${CMAKE_CURRENT_SOURCE_DIR}/libpredict.symbols")
endif()
set_target_properties(predict PROPERTIES LINK_FLAGS "${predict_link_flags}")

install(TARGETS predict predict_static
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/predict
)

# pkg-config file
configure_file(predict.pc.in predict.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/predict.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
