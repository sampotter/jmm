cmake_minimum_required (VERSION 3.15.3)

project (sjs_eikonal C)

add_library (sjs_eikonal heap.c sjs.c hermite.c vec.c)

add_executable (main main.c)
target_link_libraries (main PRIVATE sjs_eikonal)

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  add_executable (test test.c)
  target_link_libraries (test PRIVATE sjs_eikonal)
endif ()
