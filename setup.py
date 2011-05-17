#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    print '(WARNING: importing distutils, not setuptools!)'
    from distutils.core import setup

setup(name='pique',
      version='0.1',
      description='A efficient peak finder for high coverage ChIP-seq experiments, written in python.',
      author='Russell Neches',
      author_email='ryneches@ucdavis.edu',
      url='',
      packages=['pique', 'pique.tests'],
      package_data={'pique.tests': ['test_IP_fwd.txt',  \
                                    'test_IP_rev.txt',  \
                                    'test_WCE_fwd.txt', \
                                    'test_WCE_rev.txt'  ]},
      license='BSD',
      test_suite = 'nose.collector'
      )
