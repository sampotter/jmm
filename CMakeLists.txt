cmake_minimum_required (VERSION 3.15.3)

include (CheckIPOSupported)
check_ipo_supported (RESULT IPO_SUPPORTED)

project (jmm)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -Werror -mavx")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror -mavx")

set (JMM_SRCS
  bicubic.c
  bucket.c
  cubic.c
  dial.c
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

add_library (jmm STATIC ${JMM_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET jmm PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

# Masha's library

set (MASHA_LIB_SRCS
  JMMversion3.c
  mesh.c
  BucketSort.c
  HeapSort.c
  QuickSort.c
  JMMupdates.c
  linear_algebra.c
  slowness_and_uexact.c
  Newton.c
  )

add_executable (JMMversion3 ${MASHA_LIB_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET JMMversion3 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

# Unit tests

find_package (Catch2 REQUIRED)
add_executable (tests
  test_main.cpp
  test_dial.cpp
  test_vec.cpp
  )
target_link_libraries (tests Catch2::Catch2)
target_link_libraries (tests jmm)
if (IPO_SUPPORTED)
  set_property (TARGET tests PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

include (CTest)
include (Catch)
catch_discover_tests (tests)
if (IPO_SUPPORTED)
  set_property (TARGET dial PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

add_executable (main main.c)
target_link_libraries (main jmm)
if (IPO_SUPPORTED)
  set_property (TARGET main PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()
