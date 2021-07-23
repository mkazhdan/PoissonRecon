import os
import numpy

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# ------------------------------------------------------------------------------
# paths to these external libraries
libs = ['png', 'z', 'jpeg']

# TODO: need to properly pass in the paths to these external libs
ext_inc_dir = ['/opt/local/include']
ext_lib_dir = ['/opt/local/lib']

# ------------------------------------------------------------------------------
# path to the project directory
proj_dir = os.path.dirname(os.path.realpath('.'))

# code to be build for this extension module
sources = [os.path.join(proj_dir, 'Src', 'PoissonReconLib.cpp'),
           'pypoissonrecon.pyx']

include_dirs = ext_inc_dir + [proj_dir]
lib_dirs = ext_lib_dir

# ------------------------------------------------------------------------------
# flags taken from original Makefile
CXX_FLAGS = "-fopenmp -std=c++14 -pthread -fPIC -Wno-deprecated -Wno-invalid-offsetof "
CXX_LINKER_FLAGS = "-fopenmp -lstdc++ -lpthread "

 # let's just assume we are building release
if True:
    CXX_FLAGS += "-Ofast -DRELEASE -funroll-loops -ffast-math "
    CXX_LINKER_FLAGS += "-Ofast "

# if we wanted to build debug
else:
    CXX_FLAGS += "-g3 "
    CXX_LINKER_FLAGS += ""

# flags to suppress further warnings
if True:
    CXX_FLAGS += "-Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-local-typedefs "
    CXX_FLAGS += "-Wno-delete-non-virtual-dtor -Wno-class-memaccess -Wno-strict-aliasing "
    CXX_FLAGS += "-Wno-sign-compare -Wno-reorder -Wno-dangling-else "
    CXX_FLAGS += "-Wno-maybe-uninitialized -Wno-format-overflow "

# ------------------------------------------------------------------------------
extn_mod = Extension('pypoissonrecon', language='c++',
                     sources = sources,
                     include_dirs = include_dirs + [numpy.get_include()],
                     extra_compile_args = CXX_FLAGS.split(),
                     extra_link_args = CXX_LINKER_FLAGS.split(),
                     libraries=libs,
                     library_dirs=lib_dirs)

# ------------------------------------------------------------------------------
setup(
    name='pypoissonrecon',
    version='13.72',
    description='Poisson Surface Reconstruction (Adaptive Multigrid Solvers)',
    author='',
    author_email='',
    url='https://github.com/mkazhdan/PoissonRecon',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [extn_mod]
)
# ------------------------------------------------------------------------------
