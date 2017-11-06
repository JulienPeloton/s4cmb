#!/usr/bin/python
"""
Script to generate normalized inputs for the software x2pure
(pure pseudo-spectrum estimator).
For > 99.9 percent of the population, this module is useless and you do not
have to use it.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import os
import numpy as np

from s4cmb.input_sky import write_healpix_cmbmap

def safe_mkdir(path, verbose=False):
    """
    Create a path and catch the race condition between path exists and mkdir.

    Parameters
    ----------
    path : string
        Name of the folder to be created.
    verbose : bool
        If True, print messages about the status.

    Examples
    ----------
    Folders are created
    >>> safe_mkdir('toto/titi/tutu', verbose=True)

    Folders aren't created because they already exist
    >>> safe_mkdir('toto/titi/tutu', verbose=True)
    Folders toto/titi/tutu already exist. Not created.

    >>> os.removedirs('toto/titi/tutu')
    """
    abspath = os.path.abspath(path)
    if not os.path.exists(abspath):
        try:
            os.makedirs(abspath)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    else:
        if verbose:
            print("Folders {} already exist. Not created.".format(path))

def qu_weight_mineig(cc, cs, ss, epsilon=0., verbose=False):
    """
    Create a weight map using the smallest eigenvalue
    of the polarization matrix.

    Parameters
    ----------
    cc : 1d array
        Projected noise weighted cosine**2 per-pixel.
    ss : 1d array
        Projected noise weighted sine**2 per-pixel.
    cs : 1d array
        Projected noise weighted cosine * sine per-pixel.
    epsilon : float, optional
        Threshold for selecting the pixels. 0 <= epsilon < 1/4.
        The higher the more selective.

    Returns
    ----------
    weights : 1d array
        Polarisation weights per-pixel.

    """
    trace = cc + ss
    trace2 = trace * trace
    det = (cc * ss - cs * cs)

    val2 = trace2 - 4 * det
    valid = (val2 > 0.0)
    val = np.zeros_like(val2)
    val[valid] = np.sqrt(val2[valid])

    weight = np.zeros_like(trace)
    lambda_minus = (trace - val) / 2

    if verbose:
        print('criterion is', epsilon, '< det < 1/4 (epsilon= 0. by default)')

    valid = (lambda_minus > (trace - np.sqrt(trace2 - 4*epsilon * trace2)) / 2)

    if verbose:
        valid3 = [x for x in valid if x is True]
        print('number of pixels kept:', len(valid3), '/', np.sum(trace > 0))
        print('Percentage cut: {:3f} %%'.format(
            (1. - float(len(valid3)) / np.sum(trace > 0)) * 100.))

    weight[valid] = lambda_minus[valid]
    return weight

def write_maps_a_la_xpure(OutputSkyMap, name_out, output_path):
    """
    Grab sky data from OutputSkyMap and write them into files readable by
    the software xpure.

    Parameters
    ----------
    OutputSkyMap : OutputSkyMap instance
        Instance of OutputSkyMap containing map data.
    name_out : string
        File name (.fits) where to store the data.
    output_path : string
        Folder where to put the data.

    """
    safe_mkdir(os.path.join(output_path, name_out))

    fits_I = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    fits_Q = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    fits_U = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))

    I, Q, U = OutputSkyMap.get_IQU()
    fits_I[OutputSkyMap.obspix] = I
    fits_Q[OutputSkyMap.obspix] = Q
    fits_U[OutputSkyMap.obspix] = U

    full_path = os.path.join(
        output_path, name_out, 'IQU_{}.fits'.format(name_out))

    write_healpix_cmbmap(full_path,
                         data=[fits_I, fits_Q, fits_U],
                         fits_IDL=False,
                         coord=None,
                         colnames=['I', 'Q', 'U'],
                         partial=False,
                         nest=False)

    del fits_I, fits_Q, fits_U

def write_weights_a_la_xpure(OutputSkyMap, name_out, output_path, epsilon,
                             HWP=False):
    """
    Grab weight and mask from OutputSkyMap and write them into files
    readable by the software xpure.

    Parameters
    ----------
    OutputSkyMap : OutputSkyMap instance
        Instance of OutputSkyMap containing map data.
    name_out : string
        File name (.fits) where to store the data.
    output_path : string
        Folder where to put the data.
    epsilon : float
        Threshold for selecting the pixels. 0 <= epsilon < 1/4.
        The higher the more selective.
    HWP : bool
        If True, use demodulation syntax for the weights (w0, w4).
        Default is False (pair difference syntax: w, cc, ss, cs)

    """
    safe_mkdir(os.path.join(output_path, name_out))

    ## Intensity
    weight = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))

    if not HWP:
        weight[OutputSkyMap.obspix] = OutputSkyMap.w
    else:
        weight[OutputSkyMap.obspix] = OutputSkyMap.w0

    full_path = os.path.join(
        output_path, name_out, 'Iw_{}.fits'.format(name_out))
    write_healpix_cmbmap(full_path,
                         data=weight,
                         fits_IDL=False,
                         coord=None,
                         colnames=['I'],
                         partial=False,
                         nest=False)

    binary = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    mask = np.where((weight > 0))[0]
    binary[mask] = 1.0
    full_path = os.path.join(
        output_path, name_out, 'Iw_{}_norm.fits'.format(name_out))
    write_healpix_cmbmap(full_path,
                         data=binary,
                         fits_IDL=False,
                         coord=None,
                         colnames=['I'],
                         partial=False,
                         nest=False)

    ## Polarisation
    weight = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    if not HWP:
        weight[OutputSkyMap.obspix] = qu_weight_mineig(
            OutputSkyMap.cc, OutputSkyMap.cs, OutputSkyMap.ss, epsilon)
    else:
        weight[OutputSkyMap.obspix] = OutputSkyMap.w4

    full_path = os.path.join(
        output_path, name_out, 'Pw_{}.fits'.format(name_out))
    write_healpix_cmbmap(full_path,
                         data=weight,
                         fits_IDL=False,
                         coord=None,
                         colnames=['P'],
                         partial=False,
                         nest=False)

    binary = np.zeros((12 * OutputSkyMap.nside * OutputSkyMap.nside))
    mask = np.where((weight > 0))[0]
    binary[mask] = 1.0
    full_path = os.path.join(
        output_path, name_out, 'Pw_{}_norm.fits'.format(name_out))
    write_healpix_cmbmap(full_path,
                         data=binary,
                         fits_IDL=False,
                         coord=None,
                         colnames=['P'],
                         partial=False,
                         nest=False)

def create_batch(batch_file, name_out, params_s4cmb, params_xpure):
    """
    Create submission file for the software xpure

    Parameters
    ----------
    batch_file : string
        Filename.
    params_s4cmb : NormaliseS4cmbParser instance
        Object with s4cmb parameter values.
    params_xpure : NormaliseXpureParser instance
        Object with xpure parameter values.
    """
    host = os.environ['NERSC_HOST']
    if host == 'edison':
        params_xpure.nproc_per_node = 24
    elif host == 'cori':
        params_xpure.nproc_per_node = 32

    with open(batch_file, 'w') as f:
        print('#!/bin/bash -l', file=f)
        print('#SBATCH -p {}'.format(params_xpure.queue), file=f)
        print('#SBATCH -N {}'.format(params_xpure.node), file=f)
        print('#SBATCH -t {}'.format(params_xpure.time), file=f)
        print('#SBATCH -J {}'.format(name_out), file=f)
        if host == 'cori':
            print('#SBATCH -C haswell', file=f)

        print(' ', file=f)
        print('######## TO BE CHANGED BY USER ########', file=f)
        print('## Name for files and folders', file=f)
        print('name={}'.format(name_out), file=f)

        print(' ', file=f)
        print('## Radius for apodization', file=f)
        print('radius={}'.format(params_xpure.radius_apodization), file=f)

        print(' ', file=f)
        print('## Parameters for XPURE', file=f)
        print('LMAX_USER={}'.format(params_xpure.lmax_user), file=f)
        print('NSIDE_USER={}'.format(params_s4cmb.nside_out), file=f)

        print(' ', file=f)
        print('## Number of simulation (put 1 if not simulation)', file=f)
        print('NSIMU_USER=1', file=f)

        print(' ', file=f)
        print('## Number of masks (put one if several maps but one common masks)', file=f)
        print('NMASK_USER=1', file=f)

        print(' ', file=f)
        print('## Number of maps (if > 1, will do the cross correlations)', file=f)
        print('NMAPS_USER=1', file=f)

        print(' ', file=f)
        print('## Beam and bins', file=f)
        print('BEAMFILE={}'.format(params_xpure.beam_file), file=f)
        print('BINTAB={}'.format(params_xpure.bin_file), file=f)

        print(' ', file=f)
        print('## XPURE running mode (0=xpol, 1=xpure, 2=hybrid)', file=f)
        print('MODE_XPURE={}'.format(params_xpure.xpure_mode), file=f)

        print(' ', file=f)
        print('## FULL (0=myapodizemask, create_mll and XPURE) or FAST (1=only XPURE) or SEMI-FAST (=2 create_mll and XPURE)', file=f)
        print('FAST={}'.format(params_xpure.fast), file=f)

        print(' ', file=f)
        print('######## END TO BE CHANGED BY USER ########', file=f)

        print(' ', file=f)
        print('#ENVIRONEMENT', file=f)
        print('SCRATCH2=/global/cscratch1/sd/peloton', file=f)
        if host == 'edison':
            print('BINDIR=${HOME}/mapmaker/xpure/xpure/trunk/build/edison', file=f)
        elif host == 'cori':
            print('BINDIR=${HOME}/mapmaker/xpure/xpure/trunk/build/cori', file=f)
        print('THEORYDIR=${SCRATCH2}/s4cmb/additional_files', file=f)
        print('OUTPUTMLL=${SCRATCH2}/s4cmb/xpure/cls/${name}', file=f)
        print('OUTPUTCL=${SCRATCH2}/s4cmb/xpure/cls/${name}', file=f)
        print('mkdir -p ${OUTPUTMLL}', file=f)
        print('mkdir -p ${OUTPUTCL}', file=f)
        print('MASKDIR=${SCRATCH2}/s4cmb/xpure/masks/${name}', file=f)
        print('MAPDIR=${SCRATCH2}/s4cmb/xpure/maps/${name}', file=f)

        print(' ', file=f)
        print('####################################################', file=f)
        print('# MASK (COMMON MASK)', file=f)
        print('####################################################', file=f)
        print('cd $MASKDIR', file=f)

        print(' ', file=f)
        print('binary_name=$(basename $(find -name "Iw*${name}_norm.fits"))', file=f)
        print('weight_name=$(basename $(find -name "Iw*${name}.fits"))', file=f)
        print('BINARY_MASK_I1=${MASKDIR}/${binary_name}', file=f)
        print('WEIGHT_I1=${MASKDIR}/${weight_name}', file=f)
        print('APODIZED_MASK_I1=${MASKDIR}/ap30_${weight_name}', file=f)

        print(' ', file=f)
        print('binary_name=$(basename $(find -name "Pw*${name}_norm.fits"))', file=f)
        print('weight_name=$(basename $(find -name "Pw*${name}.fits"))', file=f)
        print('BINARY_MASK_P1=${MASKDIR}/${binary_name}', file=f)
        print('WEIGHT_P1=${MASKDIR}/${weight_name}', file=f)
        print('APODIZED_MASK_P1=${MASKDIR}/ap30_${weight_name}', file=f)

        print(' ', file=f)
        print('### Intensity', file=f)
        print('output_spin0_I1=${MASKDIR}/spin0_I_${name}.fits', file=f)
        print('output_spin1_I1=${MASKDIR}/spin1_I_${name}.fits', file=f)
        print('output_spin2_I1=${MASKDIR}/spin2_I_${name}.fits', file=f)

        print(' ', file=f)
        print('### Polarization', file=f)
        print('output_spin0_P1=${MASKDIR}/spin0_P_${name}.fits', file=f)
        print('output_spin1_P1=${MASKDIR}/spin1_P_${name}.fits', file=f)
        print('output_spin2_P1=${MASKDIR}/spin2_P_${name}.fits', file=f)

        print(' ', file=f)
        print('####################################################', file=f)
        print('# MAPS', file=f)
        print('####################################################', file=f)
        print('cd ${MAPDIR}', file=f)

        print(' ', file=f)
        print('for (( i=1; i<=${NMAPS_USER}; i++)); do', file=f)
        print('	MAPS[${i}]=$(basename $(find -name "IQU*${name}.fits"))', file=f)
        print('done', file=f)

        print(' ', file=f)
        print('## CHECK', file=f)
        print('for (( i=1; i<=${NMAPS_USER}; i++)); do', file=f)
        print('	echo ${MAPS[${i}]}', file=f)
        print('done', file=f)

        print(' ', file=f)
        print('cd $PBS_O_WORKDIR', file=f)

        print(' ', file=f)
        print('#########################################################################################', file=f)
        print('# From binary mask to binary apodized mask', file=f)
        print('#########################################################################################', file=f)
        print('if [ "${FAST}" -eq "0" ]', file=f)
        print('	then', file=f)
        print('	time srun -N 1 -n {} ${{BINDIR}}/myapodizemask ${{BINARY_MASK_I1}} ${{APODIZED_MASK_I1}} -minpix 1 -inside 1 -radius ${{radius}} & time srun -N 1 -n {} ${{BINDIR}}/myapodizemask ${{BINARY_MASK_P1}} ${{APODIZED_MASK_P1}} -minpix 1 -inside 1 -radius ${{radius}}'.format(params_xpure.nproc_apo, params_xpure.nproc_apo), file=f)
        print('	wait', file=f)
        print('else', file=f)
        print('	echo "Go fast - skip myapodizemask"', file=f)
        print('fi', file=f)

        print(' ', file=f)
        print('#########################################################################################', file=f)
        print('# Computation of spin window function', file=f)
        print('#########################################################################################', file=f)

        print(' ', file=f)
        print('cat > param_all_I11${name}.par << EOF', file=f)
        print('##healpix parameters', file=f)
        print('###################', file=f)
        print('nside = ${NSIDE_USER}', file=f)
        print('lmax = ${LMAX_USER} #should not exceed 2*Nside', file=f)

        print(' ', file=f)
        print('##input file parameters', file=f)
        print('#######################', file=f)
        print('maskBinary = ${BINARY_MASK_I1}', file=f)
        print('window_spin0 = ${APODIZED_MASK_I1}', file=f)
        print('inverseNoise = ${WEIGHT_I1}', file=f)

        print(' ', file=f)
        print('##output file parameters', file=f)
        print('########################', file=f)
        print('output_spin0 = ${output_spin0_I1}', file=f)
        print('output_spin1 = ${output_spin1_I1}', file=f)
        print('output_spin2 = ${output_spin2_I1}', file=f)
        print('EOF', file=f)

        print(' ', file=f)
        print('cat > param_all_P11${name}.par << EOF', file=f)
        print('##healpix parameters', file=f)
        print('###################', file=f)
        print('nside = ${NSIDE_USER}', file=f)
        print('lmax = ${LMAX_USER} #should not exceed 2*Nside', file=f)

        print(' ', file=f)
        print('##input file parameters', file=f)
        print('#######################', file=f)
        print('maskBinary = ${BINARY_MASK_P1}', file=f)
        print('window_spin0 = ${APODIZED_MASK_P1}', file=f)
        print('inverseNoise = ${WEIGHT_P1}', file=f)

        print(' ', file=f)
        print('##output file parameters', file=f)
        print('########################', file=f)
        print('output_spin0 = ${output_spin0_P1}', file=f)
        print('output_spin1 = ${output_spin1_P1}', file=f)
        print('output_spin2 = ${output_spin2_P1}', file=f)
        print('EOF', file=f)

        print(' ', file=f)
        print('if [ "${FAST}" -eq "0" ]', file=f)
        print('	then', file=f)
        print('	time srun -N {} -n {} ${{BINDIR}}/scalar2spin param_all_I11${{name}}.par >& output_scalar2spinI11${{name}} & time srun -N {} -n {} ${{BINDIR}}/scalar2spin param_all_P11${{name}}.par >& output_scalar2spinP11${{name}}'.format(int(params_xpure.nproc_scalar_to_spin // params_xpure.nproc_per_node), params_xpure.nproc_scalar_to_spin, int(params_xpure.nproc_scalar_to_spin // params_xpure.nproc_per_node), params_xpure.nproc_scalar_to_spin), file=f)
        print('	wait', file=f)
        print('else', file=f)
        print('	echo "Go fast - skip scalar2spin"', file=f)
        print('fi', file=f)

        print(' ', file=f)
        print('#########################################################################################', file=f)
        print('# X2pure', file=f)
        print('#########################################################################################', file=f)
        print('cd ${OUTPUTMLL}', file=f)

        print(' ', file=f)
        print('#########################################################################################', file=f)
        print('# Mode-mode mixing matrix', file=f)
        print('#########################################################################################', file=f)
        print('cat > createMll.par << EOF', file=f)

        print(' ', file=f)
        print('######### MODE #############', file=f)
        print('# 0 : Standard formalism', file=f)
        print('# 1 : Pure formalism', file=f)
        print('# 2 : Hybrid formalism', file=f)
        print('############################', file=f)
        print('mode = ${MODE_XPURE}', file=f)

        print(' ', file=f)
        print('############ SETUP #########', file=f)
        print('nside = ${NSIDE_USER}', file=f)
        print('lmax = ${LMAX_USER}', file=f)
        print('nmask = ${NMASK_USER}', file=f)

        print(' ', file=f)
        print('EOF', file=f)

        print(' ', file=f)
        print('for (( i=1; i<=${NMASK_USER}; i++)); do', file=f)
        print('	ind=$(($i - 1))', file=f)
        print('	cat >> createMll.par << EOF', file=f)

        print(' ', file=f)
        print('maskfile${i}_T  = ${output_spin0_I1}', file=f)

        print(' ', file=f)
        print('maskfile${i}_E_spin0 = ${output_spin0_P1}', file=f)
        print('maskfile${i}_E_spin1 = ${output_spin1_P1}', file=f)
        print('maskfile${i}_E_spin2 = ${output_spin2_P1}', file=f)

        print(' ', file=f)
        print('maskfile${i}_B_spin0 = ${output_spin0_P1}', file=f)
        print('maskfile${i}_B_spin1 = ${output_spin1_P1}', file=f)
        print('maskfile${i}_B_spin2 = ${output_spin2_P1}', file=f)

        print(' ', file=f)
        print('mllfile_TT_TT_${i} = ${OUTPUTMLL}/mll_TT_TT_BinMask${i}.fits', file=f)

        print(' ', file=f)
        print('mllfile_EE_EE_${i} = ${OUTPUTMLL}/mll_spinEE_EE_pcg${i}.fits', file=f)
        print('mllfile_EE_BB_${i} = ${OUTPUTMLL}/mll_spinEE_BB_pcg${i}.fits', file=f)
        print('mllfile_EE_EB_${i} = ${OUTPUTMLL}/mll_spinEE_EB_pcg${i}.fits', file=f)
        print('mllfile_BB_BB_${i} = ${OUTPUTMLL}/mll_spinBB_BB_pcg${i}.fits', file=f)
        print('mllfile_BB_EE_${i} = ${OUTPUTMLL}/mll_spinBB_EE_pcg${i}.fits', file=f)
        print('mllfile_BB_EB_${i} = ${OUTPUTMLL}/mll_spinBB_EB_pcg${i}.fits', file=f)

        print(' ', file=f)
        print('mllfile_TE_TE_${i} = ${OUTPUTMLL}/mll_spinTE_TE_pcg${i}.fits', file=f)
        print('mllfile_TE_TB_${i} = ${OUTPUTMLL}/mll_spinTE_TB_pcg${i}.fits', file=f)
        print('mllfile_TB_TE_${i} = ${OUTPUTMLL}/mll_spinTB_TE_pcg${i}.fits', file=f)
        print('mllfile_TB_TB_${i} = ${OUTPUTMLL}/mll_spinTB_TB_pcg${i}.fits', file=f)

        print(' ', file=f)
        print('mllfile_EB_EB_${i} = ${OUTPUTMLL}/mll_spinEB_EB_pcg${i}.fits', file=f)
        print('mllfile_EB_EE_${i} = ${OUTPUTMLL}/mll_spinEB_EE_pcg${i}.fits', file=f)
        print('mllfile_EB_BB_${i} = ${OUTPUTMLL}/mll_spinEB_BB_pcg${i}.fits', file=f)

        print(' ', file=f)
        print('EOF', file=f)

        print(' ', file=f)
        print('done', file=f)

        print(' ', file=f)
        print('if [ "${FAST}" -eq "0" ]', file=f)
        print('        then', file=f)
        print('	time srun -N {} -n {} ${{BINDIR}}/x2pure_create_mll createMll.par'.format(int(params_xpure.nproc_mll//params_xpure.nproc_per_node), params_xpure.nproc_mll), file=f)

        print('	rm -f createMll.par', file=f)
        print('elif [ "${FAST}" -eq "2" ]', file=f)
        print('        then', file=f)
        print('        time srun -N {} -n {} ${{BINDIR}}/x2pure_create_mll createMll.par'.format(int(params_xpure.nproc_mll//params_xpure.nproc_per_node), params_xpure.nproc_mll), file=f)
        print('        rm -f createMll.par', file=f)
        print('else', file=f)
        print('	echo "Go fast - skip x2pure_create_mll"', file=f)
        print('fi', file=f)

        print(' ', file=f)
        print('#########################################################################################', file=f)
        print('# X2pure', file=f)
        print('#########################################################################################', file=f)
        print('for (( n=0; n<${NSIMU_USER}; n++ )); do', file=f)

        print(' ', file=f)
        print('    echo "************************ simu $n ************************"', file=f)

        print(' ', file=f)
        print('    num=$n', file=f)

        print(' ', file=f)
        print('    #CREATE PARAMETER FILE', file=f)
        print('    cat > xpure.par << _EOF_', file=f)

        print(' ', file=f)
        print('mode = ${MODE_XPURE}', file=f)

        print(' ', file=f)
        print('nside = ${NSIDE_USER}', file=f)
        print('nmaps = ${NMAPS_USER}', file=f)
        print('nmasks = ${NMASK_USER}', file=f)

        print(' ', file=f)
        print('_EOF_', file=f)

        print(' ', file=f)
        print('    for ((j=1; j<=${NMAPS_USER}; j++)); do', file=f)

        print(' ', file=f)
        print('	cat >> xpure.par << _EOF_', file=f)

        print(' ', file=f)
        print('bellfile${j} = ${THEORYDIR}/${BEAMFILE}', file=f)
        print('mapfile${j} = ${MAPDIR}/${MAPS[${j}]}', file=f)

        print(' ', file=f)
        print('_EOF_', file=f)

        print(' ', file=f)
        print('    done', file=f)

        print(' ', file=f)
        print('    cat >> xpure.par << _EOF_', file=f)

        print(' ', file=f)
        print('lmaxSim = ${LMAX_USER}', file=f)

        print(' ', file=f)
        print('_EOF_', file=f)

        print(' ', file=f)
        print('    for(( i=1; i<=${NMASK_USER}; i++)); do', file=f)
        print('	ind=$(($i - 1))', file=f)
        print('	cat >> xpure.par << EOF', file=f)

        print(' ', file=f)
        print('maskfile${i}_T  = ${output_spin0_I1}', file=f)

        print(' ', file=f)
        print('maskfile${i}_E_spin0 = ${output_spin0_P1}', file=f)
        print('maskfile${i}_E_spin1 = ${output_spin1_P1}', file=f)
        print('maskfile${i}_E_spin2 = ${output_spin2_P1}', file=f)

        print(' ', file=f)
        print('maskfile${i}_B_spin0 = ${output_spin0_P1}', file=f)
        print('maskfile${i}_B_spin1 = ${output_spin1_P1}', file=f)
        print('maskfile${i}_B_spin2 = ${output_spin2_P1}', file=f)

        print(' ', file=f)
        print('mllfile_TT_TT_${i} = ${OUTPUTMLL}/mll_TT_TT_BinMask${i}.fits', file=f)

        print(' ', file=f)
        print('mllfile_EE_EE_${i} = ${OUTPUTMLL}/mll_spinEE_EE_pcg${i}.fits', file=f)
        print('mllfile_EE_BB_${i} = ${OUTPUTMLL}/mll_spinEE_BB_pcg${i}.fits', file=f)
        print('mllfile_EE_EB_${i} = ${OUTPUTMLL}/mll_spinEE_EB_pcg${i}.fits', file=f)
        print('mllfile_BB_BB_${i} = ${OUTPUTMLL}/mll_spinBB_BB_pcg${i}.fits', file=f)
        print('mllfile_BB_EE_${i} = ${OUTPUTMLL}/mll_spinBB_EE_pcg${i}.fits', file=f)
        print('mllfile_BB_EB_${i} = ${OUTPUTMLL}/mll_spinBB_EB_pcg${i}.fits', file=f)

        print(' ', file=f)
        print('mllfile_TE_TE_${i} = ${OUTPUTMLL}/mll_spinTE_TE_pcg${i}.fits', file=f)
        print('mllfile_TE_TB_${i} = ${OUTPUTMLL}/mll_spinTE_TB_pcg${i}.fits', file=f)
        print('mllfile_TB_TE_${i} = ${OUTPUTMLL}/mll_spinTB_TE_pcg${i}.fits', file=f)
        print('mllfile_TB_TB_${i} = ${OUTPUTMLL}/mll_spinTB_TB_pcg${i}.fits', file=f)

        print(' ', file=f)
        print('mllfile_EB_EB_${i} = ${OUTPUTMLL}/mll_spinEB_EB_pcg${i}.fits', file=f)
        print('mllfile_EB_EE_${i} = ${OUTPUTMLL}/mll_spinEB_EE_pcg${i}.fits', file=f)
        print('mllfile_EB_BB_${i} = ${OUTPUTMLL}/mll_spinEB_BB_pcg${i}.fits', file=f)

        print(' ', file=f)
        print('EOF', file=f)
        print('    done', file=f)

        print(' ', file=f)
        print('    cat >> xpure.par << _EOF_', file=f)
        print('noise_biasT_1 = 0.', file=f)
        print('noise_biasT_2_2 = 0.', file=f)
        print('noise_biasT_1_2 = 0.', file=f)

        print(' ', file=f)
        print('noise_biasP_1 = 0.', file=f)
        print('noise_biasP_2_2 = 0.', file=f)
        print('noise_biasP_1_2 = 0.', file=f)

        print(' ', file=f)
        print('bintab = ${THEORYDIR}/${BINTAB}', file=f)
        print('#mask_list = ${INPUT}/bin2mask_43bins.fits', file=f)
        print('_EOF_', file=f)

        print(' ', file=f)
        print('cat >> xpure.par << _EOF_', file=f)
        print('pseudofile = ${OUTPUTCL}/pseudopure_pcg_${name}_$num', file=f)
        print('cellfile = ${OUTPUTCL}/cellpure_pbear_${name}_$num', file=f)
        print('lmax = ${LMAX_USER}', file=f)
        print('_EOF_', file=f)
        print('    #RUN', file=f)
        print('    time srun -N {} -n {} ${{BINDIR}}/x2pure xpure.par'.format(int(params_xpure.nproc_xpure//params_xpure.nproc_per_node), params_xpure.nproc_xpure), file=f)

        print(' ', file=f)
        print('done', file=f)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
