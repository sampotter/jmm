project(
  'jmm',
  'c', 'cpp', 'cython',
  license : 'BSD-3-Clause',
  meson_version : '>=0.64.0',
  version : '0.1.1',
  default_options : [
    'c_std=c99',
    'cpp_std=c++11',
    'warning_level=2',
  ]
)

fs = import('fs')

m_dep = meson.get_compiler('c').find_library('m', required : false)
argp_dep = meson.get_compiler('c').find_library('argp', required : false)
gsl_dep = dependency('gsl')
openmp_dep = dependency('openmp')
tetgen_dep = dependency('tetgen')

jmm_lib_src = [
  'src/alist.c',
  'src/array.c',
  'src/bb.c',
  'src/bicubic.c',
  'src/bmesh.c',
  'src/bucket.c',
  'src/camera.c',
  'src/cubic.c',
  # src/dial.c
  'src/edge.c',
  'src/eik.c',
  'src/eik2g1.c',
  'src/eik2m1.c',
  'src/eik2mp.c',
  'src/eik_F3.c',
  'src/eik_F4.c',
  'src/eik_S4.c',
  'src/eik3.c',
  'src/eik3hh.c',
  'src/eik3hh_branch.c',
  'src/eik3_transport.c',
  'src/error.c',
  'src/field.c',
  'src/geom.c',
  'src/grid2.c',
  'src/grid3.c',
  'src/heap.c',
  'src/hybrid.c',
  'src/index.c',
  'src/jet.c',
  'src/log.c',
  'src/mat.c',
  'src/mesh1.c',
  'src/mesh2.c',
  'src/mesh22.c',
  'src/mesh3.c',
  'src/mesh3.cpp',
  'src/mesh_util.c',
  'src/opt.c',
  'src/par.c',
  'src/pool.c',
  'src/rtree.c',
  'src/slerp.c',
  'src/solve_cubic.c',
  'src/stats.c',
  'src/triBoxOverlap.c',
  'src/uline.c',
  'src/utetra.c',
  'src/utetra_cache.c',
  'src/utd.c',
  'src/util.c',
  'src/utri.c',
  'src/utri21.c',
  'src/utri_cache.c',
  'src/vec.c',
  'src/xfer.c'
]

jmm_inc = include_directories('include')

jmm_lib = library(
  'jmm',
  jmm_lib_src,
  dependencies : [m_dep, gsl_dep, tetgen_dep],
  include_directories : jmm_inc
)

jmm_dep = declare_dependency(link_with : jmm_lib, include_directories : jmm_inc)

subdir('examples')
subdir('wrappers')
