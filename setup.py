# -*- coding: utf-8 -*-
# Copyright 2017 Julien Peloton
# Licensed under the GPL-3.0 License, see LICENSE file for details.
from setuptools import setup, find_packages
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('s4cmb', parent_package, top_path,
                           namespace_packages=['instrument'])
    config.add_extension('scanning_strategy_f',
                         sources=['instrument/scanning_strategy_f.f90'],
                         libraries=[], f2py_options=[],
                         extra_f90_compile_args=[
                             '-ffixed-line-length-1000',
                             '-O3'],
                         extra_compile_args=[''], extra_link_args=[''],)
    return config


reqs = open('requirements.txt', 'r').read().strip().splitlines()

setup(
    name='s4cmb',
    version='0.2.0',
    url='https://github.com/JulienPeloton/s4cmb',
    download_url='',
    license='GPL-3.0',
    author='Julien Peloton',
    author_email='j.peloton@sussex.ac.uk',
    description='Simulate systematic effects in the context of CMB',
    long_description=open('README.rst', 'r').read(),
    platforms='any',
    configuration=configuration,
    packages=find_packages(),
    install_requires=reqs,
    classifiers=[
        'Development Status :: 0 - Development',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: GPL-3.0 License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Instrumentation',
        'Topic :: Simulation',
        'Topic :: CMB'
    ]
)
