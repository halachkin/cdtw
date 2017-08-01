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
    sourcefiles = ['ced/cyed.pyx']
else:
    sourcefiles = ['ced/cyed.c']

ext_modules = [
    Extension("ced.cyed",
              sourcefiles,
              include_dirs=[numpy.get_include(), 'ced/'],
              extra_compile_args=[
                  "-std=c11",
                  "-Wno-incompatible-pointer-types",
                  "-Wno-unused-variable",
                  "-Wno-absolute-value",
                  "-Wno-visibility",
                  "-Wno-#warnings"]
              )
]

setup(
    name='ced',
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    packages=find_packages(),
    version='0.0.3',
    author=['Aliaksei Halachkin', 'Yukio Fukuzawa'],
    author_email=['aliaksei.h(you know what)gmail.com', 'y (...) io * ac * nz'],
    description='Edit Distance',
    long_description=open('README.md').read(),
    url='https://github.com/yfukuzaw/ced',
    keywords=['dtw', 'dynamic time warping', 'ed', 'edit distance', 'edr', 'edit distance on real sequences'],
    install_requires=['numpy'],

)
