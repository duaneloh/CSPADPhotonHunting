from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as N

ext_modules = [Extension("geo", ["geo.pyx"]), Extension("stat", ["stat.pyx"])]
setup(name = 'myFunctions', cmdclass = {'build_ext': build_ext}, include_dirs = [N.get_include()], ext_modules = ext_modules)

#To install in place: 
#python setup.py build_ext --inplace
