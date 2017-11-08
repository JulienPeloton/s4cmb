## This examples describes how to use some basic functionalities
## of s4cmb such as:
## * simulate an input sky
## * simulate an instrument
## * simulate a scanning strategy
## * simulate TOD from the 3 previous guys
## * project back your TOD to sky maps.
from __future__ import division, absolute_import, print_function

## Initialise MPI
from mpi4py import MPI

## Import modules and routines
from s4cmb.input_sky import HealpixFitsMap

from s4cmb.instrument import Hardware

from s4cmb.scanning_strategy import ScanningStrategy

from s4cmb.tod import TimeOrderedDataPairDiff
from s4cmb.tod import OutputSkyMap
from s4cmb.tod import partial2full

from s4cmb.config_s4cmb import import_string_as_module

from s4cmb.xpure import write_maps_a_la_xpure
from s4cmb.xpure import write_weights_a_la_xpure

## Other packages needed
import os
import healpy as hp
import numpy as np
import argparse

try:
    from tqdm import *
except ImportError:
    def tqdm(x):
        """
        Do nothing. Just return x.
        """
        return x

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
        help='Tag name to identify your run. E.g. run_0.')

    ## You can also pass any new arguments, or even overwrite those
    ## from the ini file.


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
                            no_quleak=params.no_quleak)

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
        d = []
        for det in tqdm(range(inst.focal_plane.nbolometer)):
            d.append(tod.map2tod(det))

        ## Project TOD to maps
        tod.tod2map(np.array(d), sky_out_tot)

    MPI.COMM_WORLD.barrier()

    ## Coaddition over all processors.
    ## Note that all processors will then have the coadded data.
    ## If you want informations at the level of each CES (or group of),
    ## use instead:
    ## final_map = OutputSkyMap(nside=nside_out, obspix=tod.obspix)
    ## final_map.coadd_MPI(sky_out_tot, MPI=MPI)
    sky_out_tot.coadd_MPI(sky_out_tot, MPI=MPI)

    ## Check that output = input
    if rank == 0:
        nhit = partial2full(sky_out_tot.nhit,
                            sky_out_tot.obspix, sky_out_tot.nside,
                            fill_with=0.0)

        sky_out = partial2full(sky_out_tot.get_I(),
                               sky_out_tot.obspix, sky_out_tot.nside,
                               fill_with=0.0)
        mask = nhit >= 1
        assert np.all(np.abs(sky_in.I[mask] - sky_out[mask]) < 1e-9), \
            ValueError("Output not equal to input!")

        sky_out = partial2full(sky_out_tot.get_QU()[0],
                               sky_out_tot.obspix, sky_out_tot.nside,
                               fill_with=0.0)
        mask = nhit >= 10
        assert np.all(np.abs(sky_in.Q[mask] - sky_out[mask]) < 1e-9), \
            ValueError("Output not equal to input!")

        sky_out = partial2full(sky_out_tot.get_QU()[1],
                               sky_out_tot.obspix, sky_out_tot.nside,
                               fill_with=0.0)
        mask = nhit >= 10
        assert np.all(np.abs(sky_in.U[mask] - sky_out[mask]) < 1e-9), \
            ValueError("Output not equal to input!")

        ## Save data on disk into fits file for later use in xpure
        name_out = '{}_{}_{}'.format(params.tag,
                                     params.name_instrument,
                                     params.name_strategy)
        write_maps_a_la_xpure(sky_out_tot, name_out=name_out,
                              output_path='xpure/maps')
        write_weights_a_la_xpure(sky_out_tot, name_out=name_out,
                                 output_path='xpure/masks',
                                 epsilon=0.08, HWP=False)

        print("All OK! Greetings from processor 0!")

    MPI.COMM_WORLD.barrier()
