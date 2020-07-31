# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension(
        'jmm',
        ['*.pyx'],
        libraries=['jmm'],
        library_dirs=['./build']
    )
]

setup(
    name='jmm',
    ext_modules=cythonize(extensions)
)
