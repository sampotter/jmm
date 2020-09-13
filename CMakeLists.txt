cmake_minimum_required (VERSION 3.18.2)

include (CheckIPOSupported)
check_ipo_supported (RESULT IPO_SUPPORTED)

project (jmm VERSION 0.1.1)

set (INSTALL_INCLUDE_DIR include)
set (INSTALL_LIB_DIR lib)

set (CMAKE_MODULE_PATH
  "${CMAKE_SOURCE_DIR}/cmake/sanitizers-cmake/cmake"
  ${CMAKE_MODULE_PATH}
  )

find_package (Sanitizers)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -Werror -mavx")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror -mavx")

set (JMM_PUBLIC_HEADERS
  def.h
  dial.h
  )

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
  vec.c
  )

add_library (jmm STATIC ${JMM_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET jmm PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()
if (UNIX)
  target_link_libraries (jmm m)
endif ()
add_sanitizers (jmm)
target_include_directories (
  jmm INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
  )

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
if (UNIX)
  target_link_libraries (JMMversion3 m)
endif ()
add_sanitizers (JMMversion3)

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
add_sanitizers (tests)

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
add_sanitizers (main)

#
# Installation and packaging
#

include (CMakePackageConfigHelpers)

configure_package_config_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/JMMConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/JMMConfig.cmake
  INSTALL_DESTINATION ${INSTALL_LIB_DIR}/JMM/cmake
  PATH_VARS INSTALL_INCLUDE_DIR
  )

write_basic_package_version_file (
  ${CMAKE_CURRENT_BINARY_DIR}/JMMConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion
  )

install (
  TARGETS jmm
  EXPORT JMMTargets
  ARCHIVE DESTINATION lib
  )

install (FILES ${JMM_PUBLIC_HEADERS} DESTINATION include/jmm)

install (EXPORT JMMTargets
  FILE JMMTargets.cmake
  DESTINATION lib/cmake/JMM
  )

install (FILES
  "${CMAKE_CURRENT_BINARY_DIR}/JMMConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/JMMConfigVersion.cmake"
  DESTINATION lib/cmake/JMM
  )

include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${jmm_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${jmm_VERSION_MINOR}")
include (CPack)
