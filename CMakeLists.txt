cmake_minimum_required (VERSION 3.15.3)

set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

project (_sjs)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -Werror")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror")

find_package (pybind11 REQUIRED)

file (GLOB SJS_SRCS *.c *.h)

pybind11_add_module (_sjs MODULE sjs.cpp ${SJS_SRCS})
target_compile_definitions (_sjs PRIVATE USING_PYBIND11=1)
set_property (TARGET _sjs PROPERTY CXX_STANDARD 17)

add_executable (scratch scratch.cpp ${SJS_SRCS})
