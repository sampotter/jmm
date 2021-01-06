cmake_minimum_required (VERSION 3.16.3)

include (CheckIPOSupported)
check_ipo_supported (RESULT IPO_SUPPORTED)

project (jmm LANGUAGES C CXX VERSION 0.1.1)

option (BUILD_TESTS "Build unit tests")

set (CMAKE_XCODE_GENERATE_SCHEME TRUE)

set (INSTALL_INCLUDE_DIR include)
set (INSTALL_LIB_DIR lib)

set (CMAKE_MODULE_PATH
  "${CMAKE_SOURCE_DIR}/cmake/sanitizers-cmake/cmake"
  ${CMAKE_MODULE_PATH}
  )

find_package (Sanitizers)

find_package (OpenMP)
if (NOT OpenMP_C_FOUND)
  message (WARNING "Failed to find OpenMP!")
endif ()

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -Werror -mavx")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror -mavx")

set (JMM_PUBLIC_HEADERS
  def.h
  dial.h
  eik3.h
  )

set (JMM_SRCS
  array.c
  bary.c
  bb.c
  bicubic.c
  bucket.c
  cubic.c
  dial.c
  eik.c
  eik_F3.c
  eik_F4.c
  eik_S4.c
  eik3.c
  field.c
  geom.c
  grid3.c
  heap.c
  hybrid.c
  index.c
  jet.c
  mat.c
  mesh3.c
  opt.c
  utetra.c
  util.c
  utri.c
  vec.c
  xfer.c
  )

add_library (jmm STATIC ${JMM_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET jmm PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()
if (UNIX)
  target_link_libraries (jmm PRIVATE m)
endif ()
if (OpenMP_C_FOUND)
  target_link_libraries (jmm PRIVATE OpenMP::OpenMP_C)
endif ()
add_sanitizers (jmm)
target_include_directories (
  jmm INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/jmm>
  )

#
# Examples
#

# basic_solve
add_executable (basic_solve examples/basic_solve.c)
target_link_libraries (basic_solve jmm)
target_include_directories (basic_solve PRIVATE jmm)
add_sanitizers (basic_solve)
set_property (
  TARGET basic_solve
  PROPERTY XCODE_SCHEME_ARGUMENTS
  "${CMAKE_CURRENT_BINARY_DIR}/sphere_verts.bin"
  "${CMAKE_CURRENT_BINARY_DIR}/sphere_cells.bin"
  "830"
  "${CMAKE_CURRENT_BINARY_DIR}/sphere_jets.bin")

# transfer_solution
add_executable (transfer_solution examples/transfer_solution.c)
target_link_libraries (transfer_solution jmm)
target_include_directories (transfer_solution PRIVATE jmm)
add_sanitizers (transfer_solution)

#
# Unit tests
#

if (BUILD_TESTS)
  find_library (CGREEN_LIBRARY cgreen HINTS /usr/local/lib)
  find_package (GSL REQUIRED)

  message ("${CGREEN_LIBRARY}")
  set (CGREEN_INCLUDE_DIRS "/usr/local/include")

  set (JMM_TESTS_SRCS
    test_array.c
    test_bb3.c
    test_bb3tri.c
    test_bb3tet.c
    test_eik3.c
    test_opt.c
    test_mat.c
    test_mesh3.c
    test_utri.c
    )
  add_library (jmm_tests SHARED ${JMM_TESTS_SRCS})
  target_link_libraries (jmm_tests PRIVATE ${CGREEN_LIBRARY})
  target_include_directories (jmm_tests PRIVATE ${CGREEN_INCLUDE_DIRS})
  target_link_libraries (jmm_tests PRIVATE jmm)
  target_link_libraries (jmm_tests PRIVATE GSL::gsl)
  add_sanitizers (jmm_tests)

  enable_testing ()
  add_test (NAME jmm COMMAND cgreen-runner $<TARGET_FILE_NAME:jmm_tests>)
endif ()

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
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${jmm_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${jmm_VERSION_MINOR}")
include (CPack)
