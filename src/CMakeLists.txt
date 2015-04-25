add_library(predict SHARED orbit.c sdp4.c sgp4.c unsorted.c vec3.c)
add_library(predict_static STATIC orbit.c sdp4.c sgp4.c unsorted.c vec3.c)
set_target_properties(predict_static PROPERTIES OUTPUT_NAME predict)

# VERSION and SOVERSION are the same, see comment in ../CMakeLists.txt
set_target_properties(predict PROPERTIES SOVERSION ${PROJECT_VERSION})
set_target_properties(predict_static PROPERTIES SOVERSION ${PROJECT_VERSION})

# Header files needed to compile an application using libpredict
set(predict_public_header)
list(APPEND predict_public_header
  observatory.h
  orbit.h
  sdp4.h
  sgp4.h
  unsorted.h
  vec3.h
)
set_target_properties(predict PROPERTIES PUBLIC_HEADER "${predict_public_header}")
# Maybe needed on Mac OS?
#set_target_properties(predict PROPERTIES FRAMEWORK TRUE)

# Tell the linker to use the symbol version script
get_property(predict_link_flags TARGET predict PROPERTY LINK_FLAGS)
set(predict_link_flags "${predict_link_flags} -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/libpredict.symver")
set_target_properties(predict PROPERTIES LINK_FLAGS ${predict_link_flags})

target_include_directories(predict PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Make libpredict.so.<major version> symlink
set(libpredict_so_major "libpredict.so.${PROJECT_VERSION_MAJOR}")
add_custom_command(
  TARGET predict POST_BUILD 
  COMMAND ${CMAKE_COMMAND} -E create_symlink "libpredict.so.${PROJECT_VERSION}" ${libpredict_so_major}
)

install(TARGETS predict predict_static
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/predict
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${libpredict_so_major}" DESTINATION lib)
