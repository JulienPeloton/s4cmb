## This examples describes how to use some basic functionalities
## of s4cmb such as:
## * simulate an input sky
## * simulate an instrument
## * simulate a scanning strategy
## * simulate TOD from the 3 previous guys
## * project back your TOD to sky maps.
## Author: j.peloton@sussex.ac.uk
from __future__ import division, absolute_import, print_function

## Initialise MPI
from mpi4py import MPI

## Import modules and routines
from s4cmb.input_sky import HealpixFitsMap

from s4cmb.instrument import Hardware

from s4cmb.scanning_strategy import ScanningStrategy

from s4cmb.tod import TimeOrderedDataPairDiff
from s4cmb.tod import OutputSkyMap

from s4cmb.systematics import waferts_add_diffbeam
from s4cmb.systematics import get_kernel_coefficients
from s4cmb.systematics import inject_beam_ellipticity

from s4cmb.config_s4cmb import import_string_as_module

from s4cmb.xpure import create_batch
from s4cmb.xpure import write_maps_a_la_xpure
from s4cmb.xpure import write_weights_a_la_xpure

## Other packages needed
import os
import healpy as hp
import numpy as np
import argparse
import commands

def addargs(parser):
    """ Parse command line arguments for s4cmb """

    ## Defaults args - load instrument, scan and sky parameters
    parser.add_argument(
        '-inifile', dest='inifile',
        required=True,
        help='Configuration file with parameter values.')
    parser.add_argument(
        '-tag', dest='tag',
        required=True,
        help='Tag name to identify your run. E.g. run_0_crosstalk.')

    ## You can also pass any new arguments, or even overwrite those
    ## from the ini file.
    parser.add_argument(
        '-type_hwp', dest='type_hwp',
        default=None,
        help='The type of HWP that you want to mount on your instrument.')
    parser.add_argument(
        '-freq_hwp', dest='freq_hwp',
        default=None, type=float,
        help='The frequency of rotation of the HWP in Hz.')
    parser.add_argument(
        '-angle_hwp', dest='angle_hwp',
        default=None, type=float,
        help='The offset of the HWP in degree.')
    parser.add_argument(
        '-input_filename', dest='input_filename',
        required=True, nargs='+',
        help='Input fits with alms.')
    parser.add_argument(
        '-array_noise_seed', dest='array_noise_seed',
        required=True, type=int,
        help='Seed to generate noise.')
    parser.add_argument(
        '-sim_number', dest='sim_number',
        required=True, type=int,
        help='Number of the sim.')
    parser.add_argument(
        '-folder_out', dest='folder_out',
        required=True,
        help='Name of the output folder.')
    parser.add_argument(
        '-nside_in', dest='nside_in',
        required=True, type=int,
        help='Name of the output folder.')
    parser.add_argument(
        '-fwhm_in', dest='fwhm_in',
        required=True, type=float,
        help='Name of the output folder.')

    ## Arguments for beam ellipticity - see s4cmb.systematics.
    parser.add_argument(
        '-mu_beamellipticity', dest='mu_beamellipticity',
        default=10., type=float,
        help='Mean for the detector beam ellipticity distribution (in %).')
    parser.add_argument(
        '-sigma_beamellipticity', dest='sigma_beamellipticity',
        default=5., type=float,
        help='Standard deviation for the ' +
        'detector beam ellipticity distribution (in %).')
    parser.add_argument(
        '-do_diffbeamellipticity', dest='do_diffbeamellipticity',
        action='store_true',
        help='If selected, inject differential beam ' +
        'ellipticity btw detectors within a pair.')
    parser.add_argument(
        '-beam_seed', dest='beam_seed',
        default=5438765, type=int,
        help='Control the random seed used to generate beam coefficients.')


if __name__ == "__main__":
    """
    Launch the pipeline!
    """
    parser = argparse.ArgumentParser(
        description='MPI version of s4cmb')
    addargs(parser)
    args = parser.parse_args(None)

    ## Import parameters from the user parameter file
    params = import_string_as_module(args.inifile)

    ## Overwrite ini file params with params pass to the App directly
    for key in args.__dict__.keys():
        new = getattr(args, key)
        if key in params.__dict__.keys():
            old = getattr(params, key)
            if new in [None, 'None'] and old is not None:
                continue
            else:
                print("Overwriting {} with new value: {} -> {}".format(
                    key, old, new))

        setattr(params, key, new)

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size

    ##################################################################
    ## START OF THE SIMULATION
    ## Flow is the following:
    ##   sky -> instrument -> scanning strategy ->
    ##      MAP2TOD -> (systematics) -> TOD2MAP
    ##################################################################
    ## Initialise our input maps
    sky_in = HealpixFitsMap(params.input_filename,
                            fwhm_in=params.fwhm_in,
                            nside_in=params.nside_in,
                            map_seed=params.map_seed,
                            do_pol=params.do_pol,
                            verbose=params.verbose,
                            no_ileak=params.no_ileak,
                            no_quleak=params.no_quleak,
                            compute_derivatives=True)

    ## Initialise our instrument
    inst = Hardware(ncrate=params.ncrate,
                    ndfmux_per_crate=params.ndfmux_per_crate,
                    nsquid_per_mux=params.nsquid_per_mux,
                    npair_per_squid=params.npair_per_squid,
                    fp_size=params.fp_size,
                    fwhm=params.fwhm,
                    beam_seed=params.beam_seed,
                    projected_fp_size=params.projected_fp_size,
                    pm_name=params.pm_name,
                    type_hwp=params.type_hwp,
                    freq_hwp=params.freq_hwp,
                    angle_hwp=params.angle_hwp,
                    verbose=params.verbose)

    ## Initialize our scanning strategy
    scan = ScanningStrategy(nces=params.nces,
                            start_date=params.start_date,
                            telescope_longitude=params.telescope_longitude,
                            telescope_latitude=params.telescope_latitude,
                            telescope_elevation=params.telescope_elevation,
                            name_strategy=params.name_strategy,
                            sampling_freq=params.sampling_freq,
                            sky_speed=params.sky_speed,
                            ut1utc_fn=params.ut1utc_fn,
                            language=params.language)
    scan.run()

    ## Let's now generate our TOD from our input sky, instrument,
    ## and scanning strategy.
    if params.verbose:
        print("Proc [{}] doing scans".format(rank), range(
            rank, scan.nces, size))

    ## Beam ellipticity
    ## Careful, seed must be the same for all CES (beam ellipticity does not
    ## change much over time).
    sig_1, sig_2, ellip_ang = inject_beam_ellipticity(
        sigma_gaussian=inst.beam_model.sig_1[0],
        mu_beamellipticity=args.mu_beamellipticity,
        sigma_beamellipticity=args.sigma_beamellipticity,
        nbolo=inst.focal_plane.nbolometer,
        do_diffbeamellipticity=args.do_diffbeamellipticity,
        seed=args.beam_seed)
    inst.beam_model.sig_1 = sig_1
    inst.beam_model.sig_2 = sig_2
    inst.beam_model.ellip_ang = ellip_ang

    ## Collect intensity derivatives
    ## Careful with the order!
    intensity_derivatives = np.array(
        [sky_in.I,
         sky_in.dIdt, sky_in.dIdp,
         sky_in.d2Idpdt, sky_in.d2Id2t, sky_in.d2Id2p])

    ## Noise seeds
    state_for_noise = np.random.RandomState(params.array_noise_seed)
    seeds_for_noise = state_for_noise.randint(0, 1e6, scan.nces)

    for pos_CES, CESnumber in enumerate(range(rank, scan.nces, size)):
        if params.verbose:
            print("Proc [{}] with seeds ".format(rank),
                  seeds_for_noise[CESnumber], seeds_for_noise)
        tod = TimeOrderedDataPairDiff(
            inst, scan, sky_in,
            CESnumber=CESnumber,
            projection=params.projection,
            nside_out=params.nside_out,
            pixel_size=params.pixel_size,
            width=params.width,
            array_noise_level=params.array_noise_level,
            array_noise_seed=seeds_for_noise[CESnumber],
            mapping_perpair=params.mapping_perpair)

        ## Initialise map containers for each processor
        if pos_CES == 0:
            sky_out_tot = OutputSkyMap(projection=tod.projection,
                                       nside=tod.nside_out,
                                       obspix=tod.obspix,
                                       npixsky=tod.npixsky,
                                       pixel_size=tod.pixel_size)

        ## Scan input map to get TODs
        for pair in tod.pair_list:
            d = np.array([
                tod.map2tod(det) for det in pair])

            ## Compute beam kernel coefficients
            K = get_kernel_coefficients(inst.beam_model, [pair])

            ## Get spurious signal from beam mismatch
            beam_orientation = tod.return_parallactic_angle(pair[0])
            waferts_add_diffbeam(
                waferts=d,
                point_matrix=tod.point_matrix,
                beam_orientation=beam_orientation,
                intensity_derivatives=intensity_derivatives,
                diffbeam_kernels=K,
                pairlist=pair,
                spins='012')

            ## Project TOD to maps
            tod.tod2map(d, sky_out_tot)

    MPI.COMM_WORLD.barrier()

    ## Coaddition over all processors.
    ## Note that all processors will then have the coadded data.
    ## If you want informations at the level of each CES (or group of),
    ## use instead:
    ## final_map = OutputSkyMap(nside=nside_out, obspix=tod.obspix)
    ## final_map.coadd_MPI(sky_out_tot, MPI=MPI)
    sky_out_tot.coadd_MPI(sky_out_tot, MPI=MPI)

    if rank == 0:
        if params.projection == 'flat':
            name_out = '{}_{}_{}'.format(params.tag,
                                         params.name_instrument,
                                         params.name_strategy)
            sky_out_tot.pickle_me(
                '{}/sim{:03d}_{}.pkl'.format(
                    args.folder_out, args.sim_number, name_out),
                shrink_maps=False, crop_maps=2**12,
                epsilon=0., verbose=False)

        elif params.projection == 'healpix':
            ## Save data on disk into fits file for later use in xpure
            name_out = 'sim{:03d}_{}_{}_{}'.format(
                args.sim_number,
                params.tag,
                params.name_instrument,
                params.name_strategy)

            write_maps_a_la_xpure(
                sky_out_tot,
                name_out=name_out,
                output_path='xpure/maps')
            write_weights_a_la_xpure(
                sky_out_tot,
                name_out=name_out,
                output_path='xpure/masks',
                epsilon=0.08, HWP=False)

            if args.inifile_xpure is not None:
                ## Import parameters from the user parameter file
                params_xpure = import_string_as_module(args.inifile_xpure)

                batch_file = 'sim{:03d}_{}_{}_{}.batch'.format(
                    args.sim_number,
                    params.tag,
                    params.name_instrument,
                    params.name_strategy)
                create_batch(batch_file, name_out, params, params_xpure)

                qsub = commands.getoutput('sbatch ' + batch_file)
                print(qsub)

    MPI.COMM_WORLD.barrier()
