# try:
# from setuptools import setup
# except ImportError:
# #This should not be the case though
# from distutils.core import setup
from distutils.core import setup    
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os.path
import sys
import numpy.core


numpy_include_dir = os.path.join(os.path.dirname(numpy.core.__file__),'include')
#fiboheap_include_dir = os.path.join('zen','util')

ext_modules = [ 
    Extension('netcontrolz.util.matching', ['netcontrolz/util/matching.pyx'], language='c++',
                include_dirs=[numpy_include_dir], libraries=['python2.7','m','util','dl'])
]

setup(
    name = 'Network Controlz',
    #version = 'X',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    packages = [
        'netcontrolz','netcontrolz.util','netcontrolz.tests'
    ],
    package_data = {
        'netcontrolz.util' : ['*.pxd'],
    },

    # # dependencies
    setup_requires = ['zen','distribute','cython>=0.14'],
    install_requires = ['numpy>=1.6.1','matplotlib>=1.0.1'],
    #
    # # testing suite
    # test_suite = 'zen.test',
    #
    # # project metadata
    author = 'Justin Ruths',
    author_email = 'jruths@utdallas.edu',
    description = 'NetControlz is a library of network control routines specifically for Python.',
    license = 'BSD',
    url = 'https://github.com/jruths/netcontrolz',
    download_url = 'https://github.com/jruths/netcontrolz'
)
