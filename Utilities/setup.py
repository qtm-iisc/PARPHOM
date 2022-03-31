import setuptools
from numpy.distutils.core import setup, Extension

ext1 = Extension(name='bz_integration', 
                 sources=['bz_integration.f90'],
                 f2py_options=['--verbose'])

ext2 = Extension(name='neighbor',
                 sources=['neighbor_list.f90'],
                 f2py_options=['--verbose'])


setup(name = 'bz_integration',
      description = 'BZ integration using various methods',
      author = 'Shinjan Mandal',
      author_email = 'mshinjan@iisc.ac.in',
      ext_modules = [ext1,ext2]
      )
