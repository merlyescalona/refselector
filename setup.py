from setuptools import setup, find_packages
def readme():
    with open('README.md') as f:
        return f.read()
long_description=''

setup(name='ngsphy-refselector',\
    version='1.0.0',\
    description='',\
    long_description=long_description,\
    url='https://github.com/merlyescalona/refselector',\
    download_url='https://github.com/merlyescalona/refselector/blob/master/dist/refselector.1.0.0.tar.gz',\
    author='Merly Escalona',\
    author_email='merlyescalona@uvigo.es',\
    license='GNU/GPL v3',\
    packages=['ngsphy-refselector'],\
    package_dir={'ngsphy-refselector': 'ngsphy-refselector'},\
    py_modules = [\

    ],\
    install_requires=[\
        'argparse',\
        'ConfigParser',\
        'datetime',\
        'logging',\
        'numpy',\
        'setuptools',\
        'scipy',\
        'sqlite3'\
    ],\
    classifiers=[\
        'Development Status :: 4 - Beta',\
        'Intended Audience :: Education',\
        'Intended Audience :: Science/Research',\
        'Intended Audience :: Developers',\
        'Topic :: Scientific/Engineering :: Bio-Informatics',\
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',\
        'Programming Language :: Python :: 2.7'\
    ],\
    keywords='biology phylogenomics next-generation sequencing coverage targeted-sequencing',\
    python_requires='~=2.7',\
    scripts=['scripts/ngsphy-refselector'],\
    entry_points={
        'console_scripts':[\
            'ngsphy-refselector = ngsphy-refselector.__main__:main'\
        ]\
    },\
    zip_safe=False\
  )
