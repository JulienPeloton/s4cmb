## Parameter file for a fake experiment.
## The format is <name_in_small_caps> = <value> <letter>
## where is used to determine the type of the value:
## S(tring), F(loat), I(nteger), B(ool), N(one).

## Run ID
name_instrument = "fake"

#####################################################################
######################### Input sky #################################
#####################################################################
## can be two things:
##  1) Input fits file containing the sky maps to scan (maps will be loaded).
##  2) CAMB lensed cl file (.dat) containing lensed power spectra with
##      order ell, TT, EE, BB, TE (maps will be created on-the-fly).
## /!\ Option 2) is more for test purposes /!\
input_filename = "s4cmb/data/test_data_set_lensedCls.dat"

## Set do_pol to False if you want to load only intensity map.
do_pol = True

## If input_filename is a CAMB lensed cl file, the code generate maps at
## a resolution nside_in and convolved with
## a beam having this FWHM_in [arcmin].
## No effect if you provide fits file with maps instead of the CAMB file.
## /!\ this is more for test purposes /!\
fwhm_in = 0.0
nside_in = 512
map_seed = 547389

# Remove either I, Q or U to remove possible leakages
no_ileak = False
no_quleak = False

# Set it to True if you are reading a map in Galactic coordinate.
# (Planck maps for example).
ext_map_gal = False

#####################################################################
######################### Instrument ################################
#####################################################################
## creates the data used to model the instrument (see s4cmb/instrument.py).
## Focal plane parameters: Number of crate plate,
## number of MUX board per crate, number of SQUID per MUX
## and number of pair of bolometers per SQUID.
ncrate = 1
ndfmux_per_crate = 1
nsquid_per_mux = 1
npair_per_squid = 4

## The size of the focal plane [cm].
fp_size = 60.

## Beam parameters of the bolometers.
## Full Width Half Maximum [arcmin].
fwhm = 3.5

## Seed used to generate angle of rotation of beam axes.
beam_seed = 58347

## Diameter of the focal plane [degree].
projected_fp_size = 3.

# Pointing model parameters of the telescope.
pm_name = "5params"

# Polarisation angle of the bolometers and HWP. Choose among CRHWP or stepped.
type_hwp = "CRHWP"

## HWP frequency [Hz] (only if CRHWP is used).
freq_hwp = 2.

## starting position of the HWP if type_hwp=CRHWP is used
## or the daily step change of the HWP if type_hwp=stepped is used.
## [deg].
angle_hwp = 0.

## Time-domain Noise level for the whole array.
## Default are noiseless simulations (i.e. = None).
## If you want to inject white noise on-the-fly, then specify a
## noise level in [u]K.sqrt(s). Careful the units has to be the same as
## the input map! Note also that it corresponds to the polarisation level.
array_noise_level = None

## Seed used to generate random numbers to simulate noise.
## From this single seed, we generate a list of seeds
## for all detectors. Has an effect only if array_noise_level is provided.
array_noise_seed = 487587

#####################################################################
######################### Scanning strategy #########################
#####################################################################
# A scanning strategy consists in defining the site of observation
# on earth for which we will make the observation, the region
# of the sky to observe, and the schedule of observations.
# (see s4cmb/scanning_strategy.py)

## Name of a pre-defined scanning strategy to define the boundaries
## of the scan: elevation, azimuth, and time. Only available for the:
## moment name_strategy = deep_patch.
name_strategy = "deep_patch"

## Number of CES to run. Should not exceed the number of CES contained in
## the defined scanning strategy set in name_strategy.
nces = 12

## Starting date for observations (The format is: YYYY/M/D HH:MM:SS)
start_date = "2013/1/1 00:00:00"

## Position of the telescope on Earth
telescope_longitude = "-67:46.816"
telescope_latitude = "-22:56.396"

## Elevation (Alt) of the telescope [meter]
telescope_elevation = 5200.

## Sampling frequency at which time samples are recorded [Hz]
sampling_freq = 8.

## Azimuth speed of the telescope [deg/s]
sky_speed = 0.4

## File containing time correction to UTC
ut1utc_fn = "s4cmb/data/ut1utc.ephem"

#####################################################################
######################### Output sky maps ###########################
#####################################################################
## Choose the projection of the output maps: healpix of flat
projection = "healpix"

## If projection=healpix, choose the nside
nside_out = 512

## If projection=flat, choose the pixel_size
pixel_size = None

## Width of the output map in degree.
## If width < len(input scan), you will get a cropped map.
## if width > len(input scan), you will
## have your entire scan + plenty of zeros.
## A choice has to be made.
width = 20.

## Project time-domain data detector-by-detector to save a lot of memory.
## Use this only if you do not have detector-to-detector correlation!
mapping_perpair = False

#####################################################################
######################### Misc ######################################
#####################################################################
## language use to perform core computations.
language = "fortran"

## Verbose mode.
verbose = False
