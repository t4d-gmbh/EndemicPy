#!/usr/bin/env python

import os
from re import MULTILINE, search
import sys

from codecs import open

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

packages = [
    'endemic',
    'endemic/nw_construct',
    'endemic/nw_spread'
]

requires = ['numpy']# 'networkx']

version = ''
with open('endemic/__init__.py', 'r') as fd:
    version = search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        fd.read(), MULTILINE
    ).group(1)

if not version:
    raise RuntimeError('Cannot find version information')

with open('README.rst', 'r', 'utf-8') as f:
    readme = f.read()
#with open('HISTORY.rst', 'r', 'utf-8') as f:
#    history = f.read()

setup(
    name='endemic',
    version=version,
    description='Simulation of disease outbreak on various host structures',
    long_description=readme + '\n\n',  # + history,
    author='Jonas I. Liechti',
    author_email='jon.liechti@gmail.com',
    url='https://tb-git.usys.ethz.ch/j-i-l/EndemicPy',
    download_url='https://tb-git.usys.ethz.ch/j-i-l/EndemicPy/tarball/%s' % str(version),
    keywords=['disease spread', 'network', 'SIS', 'SIR'],
    packages=packages,
    package_data={'': ['LICENSE', 'HISTORY']},
    package_dir={'endemic': 'endemic'},
    include_package_data=True,
    install_requires=requires,
    license='Apache 2.0',
    zip_safe=False,
    classifiers=(
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Operating System :: OS Independent'
    ),
    # dependency_links = [
    #     "http://..."
    # ],
)
