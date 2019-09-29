cmake_minimum_required (VERSION 3.15.3)

project (sjs_eikonal C)

find_library (BLOCKS_RUNTIME_LIBRARY NAMES BlocksRuntime)

add_library (sjs_eikonal heap.c sjs.c hermite.c vec.c)
target_link_libraries (sjs_eikonal ${BLOCKS_RUNTIME_LIBRARY})

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fblocks -lm")

add_executable (main main.c)
target_link_libraries (main PRIVATE sjs_eikonal)

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  add_executable (test test.c)
  target_link_libraries (test PRIVATE sjs_eikonal ${BLOCKS_RUNTIME_LIBRARY})
endif ()
