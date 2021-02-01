# -*- coding: utf-8 -*-
# Copyright 2016-2021 Julien Peloton, Giulio Fabbian
# Licensed under the GPL-3.0 License, see LICENSE file for details.
from s4cmb import __version__
from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('s4cmb', parent_package, top_path)
    config.add_extension(
        'scanning_strategy_f',
        sources=['s4cmb/scanning_strategy_f.f90'],
        libraries=[], f2py_options=[],
        extra_f90_compile_args=[
            '-ffixed-line-length-1000',
            '-O3'
        ],
        extra_compile_args=[''],
        extra_link_args=['']
    )
    config.add_extension(
        'detector_pointing_f',
        sources=['s4cmb/detector_pointing_f.f90'],
        libraries=[], f2py_options=[],
        extra_f90_compile_args=[
            '-ffixed-line-length-1000',
            '-O3'
        ],
        extra_compile_args=[''],
        extra_link_args=[''],
    )
    config.add_extension(
        'tod_f',
        sources=['s4cmb/tod_f.f90'],
        libraries=[], f2py_options=[],
        extra_f90_compile_args=[
            '-ffixed-line-length-1000',
            '-O3'
        ],
        extra_compile_args=[''],
        extra_link_args=[''],
    )
    config.add_extension(
        'systematics_f',
        sources=['s4cmb/systematics_f.f90'],
        libraries=[], f2py_options=[],
        extra_f90_compile_args=[
            '-ffixed-line-length-1000',
            '-O3'],
        extra_compile_args=[''],
        extra_link_args=[''],
    )
    return config


if __name__ == "__main__":
    # Download url
    d_url = 'https://github.com/JulienPeloton/s4cmb/archive/{}.tar.gz'.format(
        __version__)

    reqs = open('requirements.txt', 'r').read().strip().splitlines()

    package_data = {
        's4cmb/data': [
            'test_data_set_lensedCls.dat',
            'ut1utc.ephem'
        ],
        'examples': [
            'simple_app.py',
            'so_app.py',
            'so_crosstalk_app.py',
            'so_gain_variation.py',
            'simple_parameters.ini',
            'so_parameters.ini',
            'xpure_parameters.ini',
            'so_flat_MC_app.py',
            'nersc_cori.batch'
        ]
    }

    setup(
        configuration=configuration,
        version=__version__,
        url='https://github.com/JulienPeloton/s4cmb',
        download_url=d_url,
        license='GPL-3.0',
        author='Julien Peloton, Giulio Fabbian',
        author_email='peloton@lal.in2p3.fr',
        description='Simulate systematic effects in the context of CMB',
        platforms='any',
        packages=find_packages(),
        include_package_data=True,
        package_data=package_data,
        install_requires=reqs,
        classifiers=[
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Astronomy'
        ]
    )
