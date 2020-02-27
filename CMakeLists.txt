cmake_minimum_required (VERSION 3.15.3)

project (sjs_eik)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -Werror")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror")

find_package (pybind11 REQUIRED)

file (GLOB SJS_EIK_SRCS *.c *.h)

message (${SJS_EIK_SRCS})

pybind11_add_module (sjs_eik MODULE sjs_eik.cpp ${SJS_EIK_SRCS})
target_compile_definitions (sjs_eik PRIVATE USING_PYBIND11=1)
set_property (TARGET sjs_eik PROPERTY CXX_STANDARD 17)
