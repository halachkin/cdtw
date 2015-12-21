"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
# To use a consistent encoding
from codecs import open
from os import path
import numpy

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

sourcefiles = ['cdtw/cydtw.pyx']
ext_modules = [
        Extension("cdtw.cydtw", 
            sourcefiles,
            include_dirs=[numpy.get_include(), 'cdtw/'],
            # extra_compile_args=[""]
            extra_compile_args=["-std=c11"]
            )
        ]
#ext_modules = [
#        Extension('cdtw.cydtw', ['cdtw/cydtw.pyx'], include_dirs=[numpy.get_include(), '.'])
#        ]
setup(
    name='cdtw',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    packages = find_packages(),
    version='0.0.1',

    description='DTW computation',
    long_description=long_description,

    url='https://github.com/johanneshk/cdtw',


    install_requires=['cython', 'numpy'],

)
