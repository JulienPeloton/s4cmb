## Inputs for the software x2pure (pure pseudo-spectrum estimator).
## For > 99.9 percent of the population, this is useless and you do not
## have to touch or use this.

#####################################################################
######################### Xpure #####################################
#####################################################################
## Walltime and queue
time = "00:10:00"
queue = "regular"

## Nodes and processors (node is the total number of nodes, the other are
## number of processors to use).
node = 6
nproc_apo = 1
nproc_scalar_to_spin = 96
nproc_mll = 192
nproc_xpure = 96

## Few parameters to estimate the spectra
## Radius of the apodization [arcmin]
radius_apodization = 30

## Maximum multipole for the reconstruction.
lmax_user = 3000

## Running mode (0=xpol, 1=xpure, 2=hybrid)
xpure_mode = 1

## FULL (0=myapodizemask, create_mll and XPURE) or
## FAST (1=only XPURE) or SEMI-FAST (=2 create_mll and XPURE)
fast = 0

## Beam and bin file (just the name, not the full path)
beam_file = "beam_freq150_large_v2.fits"
bin_file = "bins_40to6143_step50_so.fits"
