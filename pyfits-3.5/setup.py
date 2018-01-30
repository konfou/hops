#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup


import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))
os.chdir('stsci.distutils-0.3.7')
os.system('python setup.py install')
os.chdir('..')

setup(
    name='pyfits',
    setup_requires=['d2to1>=0.2.12', 'stsci.distutils>=0.3'],
    d2to1=True,
    zip_safe=False
)
