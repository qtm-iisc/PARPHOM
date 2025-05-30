import setuptools

# Read README safely
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='pyphutil',
    version='1.0',
    description='Utility functions for PARPHOM',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Shinjan Mandal',
    author_email='shinjanm@umich.edu',
    url='https://github.com/qtm-iisc/PARPHOM',
    packages=['pyphutil', 'pyphutil.fortran_routines'],
    package_dir={'pyphutil': 'pyphutil'},
    include_package_data=True,
    zip_safe=False,
)
# Note: Fortran extension build is handled by f2py before install.
