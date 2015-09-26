link_directories(${LIBPREDICT_LIBRARY_DIR})

add_executable(dummy-t dummy-t.c)
add_test(NAME dummy COMMAND dummy-t)

add_executable(link-t link-t.c)
target_link_libraries(link-t predict)
add_test(NAME link COMMAND link-t)
