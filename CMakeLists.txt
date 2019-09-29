cmake_minimum_required (VERSION 3.15.3)

project (sjs_eikonal C)

if (APPLE)
    # Fix linking on macOS 10.14+. See:
    # https://stackoverflow.com/questions/54068035
    link_directories(/usr/local/lib)
endif ()

add_library (sjs_eikonal heap.c sjs.c hermite.c vec.c)

add_executable (main main.c)
target_link_libraries (main PRIVATE sjs_eikonal)

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  find_package (cgreen)
  add_library (test SHARED test.c)
  target_include_directories (test PRIVATE ${CGREEN_INCLUDE_DIRS})
  target_link_libraries (test PRIVATE sjs_eikonal ${CGREEN_LIBRARIES})
endif ()
