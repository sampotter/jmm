cmake_minimum_required (VERSION 3.15.3)

include (CheckIPOSupported)
check_ipo_supported (RESULT IPO_SUPPORTED)

project (_sjs)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -Werror")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror")

find_package (pybind11 REQUIRED)

set (SJS_SRCS
  bicubic.c
  cubic.c
  eik.c
  eik_F3.c
  eik_F4.c
  eik_S4.c
  field.c
  heap.c
  hybrid.c
  index.c
  mat.c
  math.c
  npy.c
  vec.c
  )

add_library (sjs STATIC ${SJS_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET sjs PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

pybind11_add_module (_sjs MODULE sjs.cpp)
target_link_libraries (_sjs PRIVATE sjs)
target_compile_definitions (_sjs PRIVATE USING_PYBIND11=1)
set_property (TARGET _sjs PROPERTY CXX_STANDARD 17)

add_executable (scratch scratch.cpp)
target_link_libraries (scratch PRIVATE sjs)

# Masha's library

set (MASHA_LIB_SRCS
  JMMversion2.c
  HeapSort.c
  QuickSort.c
  JMMupdates.c
  linear_algebra.c
  slowness_and_uexact.c
  Newton.c
  )

add_executable (JMMversion2 ${MASHA_LIB_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET JMMversion2 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()
