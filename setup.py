# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import glob
import versioneer

from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension(
        'jmm',
        ['src/jmm.pyx'],
        libraries=['jmm'],
        library_dirs=['./build']
    )
]

setup(
    name='jmm',
    version=versioneer.get_version(),
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            'embedsignature': True,
            'language_level': 3
        }
    ),
    package_dir={'': 'src'},
    cmdclass=versioneer.get_cmdclass()
)
