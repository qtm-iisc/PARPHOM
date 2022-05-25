#import setuptools
from numpy.distutils.core import setup, Extension

ext1 = Extension(name='bz_integration', 
                 sources=['pyphutil/fortran_routines/bz_integration.f90'],
                 )

ext2 = Extension(name='neighbor',
                 sources=['pyphutil/fortran_routines/neighbor_list.f90'],
                 )


setup(name = 'pyphutil',
      version = "1.0",
      description = 'Utility frunctions for Moire Phonons',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      author = 'Shinjan Mandal',
      author_email = 'mshinjan@iisc.ac.in',
      url = "https://github.com/ShinjanM/MoirePhonons",
      packages = ['pyphutil'],
      package_dir = {'pyphutil':'pyphutil/'},
      #install_requires = open('requirements.txt').read(),
      ext_modules = [ext1,ext2]
      )
