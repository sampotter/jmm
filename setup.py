# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import multiprocessing

from setuptools import Extension, find_packages, setup
from Cython.Build import cythonize

extensions = [
    Extension(
        '*',
        ['./src/jmm/*.pyx'],
        extra_compile_args=['-std=c99', '-O0'],
        include_dirs=['./src'],
        libraries=['jmm', 'gsl'],
        library_dirs=['./build'],
        language='c'
    )
]

setup(
    name='jmm',
    version='0.1.1',
    packages=['jmm'],
    package_dir={'': 'src'},
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            'embedsignature': True,
            'language_level': 3
        },
        nthreads=multiprocessing.cpu_count()
    ),
    zip_safe=False
)
