from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

sourcefiles = ['cydtw.pyx']
ext_modules = [Extension("cydtw", 
                          sourcefiles,
                          include_dirs=[numpy.get_include()],
                          # extra_compile_args=[""]
			 			  extra_compile_args=["-std=c11"]
                          )]

setup(
  name = 'cdtw app',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)