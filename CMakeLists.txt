cmake_minimum_required (VERSION 3.15.3)

project (sjs_eikonal C)

add_executable (sjs_eikonal main.c heap.c sjs.c hermite.c vec.c)
