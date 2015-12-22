from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
from codecs import open
import numpy

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

if use_cython:
    sourcefiles = ['cdtw/cydtw.pyx']
else:
    sourcefiles = ['cdtw/cydtw.c']

ext_modules = [
        Extension("cdtw.cydtw", 
            sourcefiles,
            include_dirs=[numpy.get_include(), 'cdtw/'],
            # extra_compile_args=[""]
            extra_compile_args=["-std=c11"]
            )
        ]

setup(
    name='cdtw',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    packages = find_packages(),
    version='0.0.1',
    author = 'Aliaksei Halachkin',
    author_email = 'aliaksei.h(you know what)gmail.com',
    description='DTW computation',
    long_description=open('README.md').read(),

    url='https://github.com/honeyext/cdtw',

    keywords = ['dtw', 'dynamic time warping'],
    install_requires=['numpy'],

)
