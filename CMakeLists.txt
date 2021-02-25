cmake_minimum_required (VERSION 3.12)

include (CheckIPOSupported)
check_ipo_supported (RESULT IPO_SUPPORTED)

project (jmm LANGUAGES C VERSION 0.1.1)

set (CMAKE_C_STANDARD 99)
set (CMAKE_C_STANDARD_REQUIRED TRUE)

option (BUILD_DOCS "Build Sphinx documentation." OFF)
option (BUILD_TESTS "Build unit tests" OFF)

set (CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/sanitizers-cmake/cmake"
  ${CMAKE_MODULE_PATH}
  )

# TODO: should do this in a more portable way... i.e. check for
# MacPorts...
include_directories (/opt/local/include)

find_package (Sanitizers)

set (CMAKE_XCODE_GENERATE_SCHEME TRUE)

set (INSTALL_INCLUDE_DIR include)
set (INSTALL_LIB_DIR lib)

set (EXTRA_FLAGS "${EXTRA_FLAGS} -fPIC -Wall -Wextra -Werror -mavx")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}")

set (JMM_PUBLIC_HEADERS
  src/def.h
  src/dial.h
  src/eik3.h
  )

set (JMM_SRCS
  src/alist.c
  src/array.c
  src/bb.c
  src/bicubic.c
  src/bmesh.c
  src/bucket.c
  src/cubic.c
  src/dial.c
  src/edge.c
  src/edgemap.c
  src/eik.c
  src/eik_F3.c
  src/eik_F4.c
  src/eik_S4.c
  src/eik3.c
  src/field.c
  src/geom.c
  src/grid3.c
  src/heap.c
  src/hybrid.c
  src/index.c
  src/jet.c
  src/log.c
  src/mat.c
  src/mesh2.c
  src/mesh3.c
  src/opt.c
  src/par.c
  src/pool.c
  src/rtree.c
  src/stats.c
  src/triBoxOverlap.c
  src/utetra.c
  src/util.c
  src/utri.c
  src/vec.c
  src/xfer.c
  )

add_library (jmm STATIC ${JMM_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET jmm PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()
if (UNIX)
  target_link_libraries (jmm PRIVATE m)
endif ()
target_include_directories (
  jmm INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/jmm>
  )

#
# Examples
#

# basic_solve
add_executable (basic_solve examples/basic_solve.c)
add_sanitizers (basic_solve)
target_link_libraries (basic_solve jmm)
target_include_directories (basic_solve PRIVATE jmm)
set_property (
  TARGET basic_solve
  PROPERTY XCODE_SCHEME_ARGUMENTS
  "${CMAKE_CURRENT_BINARY_DIR}/sphere_verts.bin"
  "${CMAKE_CURRENT_BINARY_DIR}/sphere_cells.bin"
  "830"
  "${CMAKE_CURRENT_BINARY_DIR}/sphere_jets.bin")

# raytracing
add_executable (raytracing examples/raytracing.c)
add_sanitizers (raytracing)
target_link_libraries (raytracing jmm)
target_include_directories (raytracing PRIVATE jmm)

# transfer_solution
add_executable (transfer_solution examples/transfer_solution.c)
add_sanitizers (transfer_solution)
target_link_libraries (transfer_solution jmm)
target_include_directories (transfer_solution PRIVATE jmm)

# viz_error_using_raytracing
add_executable (viz_error_using_raytracing examples/viz_error_using_raytracing.c)
add_sanitizers (viz_error_using_raytracing)
target_link_libraries (viz_error_using_raytracing jmm)
target_include_directories (viz_error_using_raytracing PRIVATE jmm)

#
# Unit tests
#

if (BUILD_TESTS)
  find_library (CGREEN_LIBRARY cgreen HINTS /usr/local/lib)
  find_package (GSL REQUIRED)

  set (CGREEN_INCLUDE_DIRS "/usr/local/include")

  set (JMM_TESTS_SRCS
    test/test_alist.c
    test/test_array.c
    test/test_bb3.c
    test/test_bb3tri.c
    test/test_bb3tet.c
    test/test_edgemap.c
    test/test_eik3.c
    test/test_geom.c
    test/test_opt.c
    test/test_mat.c
    test/test_mesh2.c
    test/test_mesh3.c
    test/test_utri.c
    test/test_vec.c
    )
  add_library (jmm_tests SHARED ${JMM_TESTS_SRCS})
  add_sanitizers (jmm_tests)
  target_link_libraries (jmm_tests PRIVATE ${CGREEN_LIBRARY})
  target_include_directories (jmm_tests PRIVATE ${CGREEN_INCLUDE_DIRS})
  target_link_libraries (jmm_tests PRIVATE jmm)
  target_include_directories (jmm_tests PRIVATE jmm)
  target_link_libraries (jmm_tests PRIVATE GSL::gsl)

  enable_testing ()
  add_test (NAME jmm COMMAND cgreen-runner $<TARGET_FILE_NAME:jmm_tests>)
endif ()

#
# Build documentation
#

if (BUILD_DOCS)
  find_package (Doxygen REQUIRED)
  find_package (Sphinx REQUIRED)

  set (DOXYGEN_INDEX_FILE ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  set (DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set (DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  set (SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
  set (SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR})

  configure_file (${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

  add_custom_command(
    OUTPUT ${DOXYGEN_INDEX_FILE}
    DEPENDS ${JMM_PUBLIC_HEADERS}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
    MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
    COMMENT "Generating Doxygen documentation"
    VERBATIM
    )

  configure_file (
    ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/conf.py
    @ONLY
    )

  add_custom_target (
    Sphinx ALL
    COMMAND
    ${SPHINX_EXECUTABLE}
    -b html -Dbreathe_projects.CatCutifier=${CMAKE_CURRENT_BINARY_DIR}
	-c ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR} # sourcedir
	${CMAKE_CURRENT_BINARY_DIR} # outputdir
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating Sphinx documentation"
    )

  add_custom_target (Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
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
