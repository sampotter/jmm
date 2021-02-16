# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import glob

from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension(
        'jmm',
        ['src/jmm.pyx'],
        libraries=['jmm', 'gc'],
        library_dirs=['./build', '/opt/local/lib']
    )
]

setup(
    name='jmm',
    version='0.1.1',
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            'embedsignature': True,
            'language_level': 3
        }
    )
)
