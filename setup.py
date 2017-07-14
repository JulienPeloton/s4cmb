# -*- coding: utf-8 -*-
# Copyright 2017 Julien Peloton
# Licensed under the GPL-3.0 License, see LICENSE file for details.

## A bit wacky... Wheel perhaps?
from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('s4cmb', parent_package, top_path)
    config.add_extension('scanning_strategy_f',
                         sources=['s4cmb/scanning_strategy_f.f90'],
                         libraries=[], f2py_options=[],
                         extra_f90_compile_args=[
                             '-ffixed-line-length-1000',
                             '-O3'],
                         extra_compile_args=[''], extra_link_args=[''],)
    config.add_extension('detector_pointing_f',
                         sources=['s4cmb/detector_pointing_f.f90'],
                         libraries=[], f2py_options=[],
                         extra_f90_compile_args=[
                             '-ffixed-line-length-1000',
                             '-O3'],
                         extra_compile_args=[''], extra_link_args=[''],)
    config.add_extension('tod_f',
                         sources=['s4cmb/tod_f.f90'],
                         libraries=[], f2py_options=[],
                         extra_f90_compile_args=[
                             '-ffixed-line-length-1000',
                             '-O3'],
                         extra_compile_args=[''], extra_link_args=[''],)
    return config


if __name__ == "__main__":
    reqs = open('requirements.txt', 'r').read().strip().splitlines()

    package_data = {'s4cmb': ['data']}

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
        package_data=package_data,
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
