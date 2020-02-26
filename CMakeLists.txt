cmake_minimum_required (VERSION 3.15.3)

project (sjs_eikonal C)

if (LINUX)
  find_library (BLOCKS_RUNTIME_LIBRARY NAMES BlocksRuntime)
endif ()

add_library (sjs_eikonal sjs_eikonal.c heap.c hermite.c vec.c)
target_compile_options (sjs_eikonal PRIVATE -lm)
if (LINUX)
  target_compile_options (sjs_eikonal PRIVATE -fblocks)
  target_link_libraries (sjs_eikonal PRIVATE ${BLOCKS_RUNTIME_LIBRARY})
endif ()

add_executable (main main.c)
target_link_libraries (main PRIVATE sjs_eikonal)

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  add_executable (test test.c)
  target_link_libraries (test PRIVATE sjs_eikonal)
  if (LINUX)
    target_compile_options (test PRIVATE -fblocks)
    target_link_libraries (test PRIVATE ${BLOCKS_RUNTIME_LIBRARY})
  endif ()
endif ()
