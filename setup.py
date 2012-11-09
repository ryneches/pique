#!/usr/bin/env python
try:
    from setuptools import setup, find_packages
    from setuptools import Extension 
    from Cython.Distutils import build_ext
except ImportError:
    print '(WARNING: importing distutils, not setuptools!)'
    from distutils.core import setup

ext_modules = [ Extension( 'pique.mapmaker', [ 'pique/mapmaker.pyx' ] ),
                Extension( 'pique.peak'    , [ 'pique/peak.pyx'     ] ) ]

setup(name = 'pique',
    version = '0.1',
    description = 'An efficient peak finder for high coverage ChIP-seq experiments.',
    #long_description=read('README'),
    author = 'Russell Neches',
    author_email = 'ryneches@ucdavis.edu',
    url = 'https://github.com/ryneches/pique',
    packages = ['pique', 'pique.tests'],
    package_data = {'pique.tests': ['test_IP_fwd.txt',  \
                                    'test_IP_rev.txt',  \
                                    'test_WCE_fwd.txt', \
                                    'test_WCE_rev.txt'  ]},
    license = 'BSD',
    test_suite = 'nose.collector',
    cmdclass = { 'build_ext' : build_ext },
    ext_modules = ext_modules,
    classifiers = [
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Cython',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: BSD License'
    ],
    requires = [ 'numpy', 'scipy', 'cython', 'pysam' ],
    scripts = [ 'scripts/piqueTk', 'scripts/pique' ]
    )
