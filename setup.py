from distutils.core import setup, Extension
from glob import glob
from Cython.Build import cythonize

ext = Extension('pysjs', ['sjs.pyx'], libraries=['sjs_eikonal'])

setup(name='pysjs', ext_modules=cythonize([ext]))
