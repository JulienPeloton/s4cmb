---
title: 'Simulating instrumental systematics of Cosmic Microwave Background experiments with s4cmb'
tags:
  - Python
  - cosmology
  - CMB
  - telescopes
  - systematics
authors:
  - name: Giulio Fabbian^[Corresponding author.]
    orcid: 0000-0002-3255-4695
    affiliation: "1, 3" # (Multiple affiliations must be quoted)
  - name: Julien Peloton^[Corresponding author.]
    orcid: 0000-0002-8560-4449
    affiliation: 2
affiliations:
 - name: Department of Physics & Astronomy, University of Sussex, Brighton BN1 9QH, UK
   index: 1
 - name: Université Paris-Saclay, CNRS/IN2P3, IJCLab, Orsay, France
   index: 2
 - name: School of Physics and Astronomy, Cardiff University, The Parade, Cardiff, CF24 3AA, UK
   index: 3
date: "`r format(Sys.time(), '%d %B %Y')`"
bibliography: s4cmb_paper.bib


---

# Summary

The observation of cosmic microwave background (CMB) anisotropies is one of the key probes of the standard cosmological model [@hu-dodelson2002].
The weak CMB polarization signal in particular can provide a new window to study the process of the growth of structures in the universe (galaxies, galaxy clusters etc.) through weak gravitational lensing as well as on the physics of the inflationary epoch in the primordial universe.
The low amplitude of CMB polarization has pushed CMB science toward the construction of increasingly sensitive experiments observing in multiple frequencies and employing telescopes with complex optical designs and focal planes with thousands of bolometric detectors operating in cryogenic environments [@Staggs_2018].
To fully utilize the sensitivity of these instruments for cosmology, the instrumental systematic effects must be well-characterized, understood and mitigated in the instrument design and through the analysis of the acquired data.


# Statement of need

The `s4cmb` (systematics for CMB) package is a Python package developed to study
the impact of instrumental systematic effects on measurements of CMB experiments based on bolometric detector technology.
`s4cmb` provides a unified framework to simulate raw data streams in the time domain (TODs) acquired by CMB experiments scanning the sky, and to inject in these
realistic instrumental systematics effect.
The development of `s4cmb` is built on experience and needs of the analysis of
data of the Polarbear ground-based experiment [@pb2014; @Ade:2017uvt].
It is designed to analyze real data, to guide the design of future instruments that require the estimation of specific systematic effects as well as to increase the realism of simulated data sets required in the development of data analysis methods.
The package has already been used in a number of scientific [@Puglisi:2018txk; @mirmelstein2020] and technical publications [@Salatino:2018voz; @Crowley:2018eib; @Gallardo:2018rix; @Bryan:2018mva].
It adopts several commonly used libraries in astronomy (`astropy` [@astropy], `healpy` [@healpy], `ephem` [@pyephem], `pyslalib` [@pyslalib]) and uses functions based on low-level languages wrapped in Python (e.g. fortran with `f2py`) for speeding up the most critical part of the code without losing the flexibility provided by a simple python user-friendly interface.
`s4cmb` is designed to be employed on systems of varying scale, from laptops to parallel supercomputing platforms thanks to its internal Message Passing Interface (MPI) support [@10.1145/169627.169855].
We also support packaging the entire application into a Docker container for portability.
The simplicity of the `s4cmb` framework allows to easily add new instrumental systematics to be simulated according to the users' needs.
As far as we know, s4cmb is the only dedicated package that enables the study of a wide range of instrumental simulations, from the instrument to the sky map, while being publicly available. For more general purposes, including some instrumental systematic effect simulations, users might also consider the use of `TOAST` [@theodore_kisner_2020_4270476], a software framework to simulate and process timestream data collected by telescopes focusing on efficient TOD manipulation on massively parallel architectures.

# Package structure and functionalities

One of the key feature of `s4cmb` is to be able to make robust simulations while being fast to run and easy to use. Simulations using ${\cal{O}}(10^3)$ detectors observing the sky for 5 hours with a data sampling rate of $100$Hz and reconstructing a sky map can be run on a single core in less than 10 minutes on a laptop. The modules implementing the major functional blocks of the library are (see \autoref{fig:s4cmb}):

![Schematic structure of `s4cmb`. Objects defined by the user are marked in red. From these, TODs are generated for the duration of each observation described by the scanning strategy. These can then be modified by introducing instrumental effects. The output of the code (blue) are the perturbed TODs or sky maps reconstructed from TODs to be employed in subsequent analysis steps (e.g. the computation of their angular power spectrum).\label{fig:s4cmb}](presentation.png){width=100%}

* `instrument.py` contains the class describing the CMB instrument in terms of position of its detectors in the focal plane, their optical beam shape, wiring in the readout electronics. It supports the most common focal plane designs employing multifrequency detectors and polarization modulation hardware such as stepped or continuously rotating half-wave plates (HWPs).
* `scanning_strategy.py` describes the schedule of the instrument observations and the motion of the instrument in terms of azimuth and elevation position at the telescope location on Earth as a function of time.
The schedule is divided in minimal units (scans) during which a telescope motion is repeated for its given duration (e.g. a back and forth motion for a fixed distance in azimuth at a constant elevation). The code parallelisation is done over scans, so that increasing proportionally the number of scans and the number of MPI workers keeps the runtime constant.
* `input_sky.py` contains the class describing the input sky signal model in HEALPix pixelization [@healpix]. This can include multiple components (CMB, Galactic emissions etc.) and can be read from an external file or be synthesized on the fly from its harmonic coefficients. The input sky needs be convolved with the main optical beam of the instrument $B$ so that the observed sky Stokes parameters $X^{\rm obs}\in\{I, Q, U\}$ are related to their true value on sky a $X^{\rm obs} \equiv B \circledast X$.
* `tod.py` includes the class to generate and handle TODs from an input sky signal map, an instrument design and a scanning strategy together with basic tools to reconstruct a sky map from TODs in order to mimic self consistently the basic data analysis pipeline of CMB experiments.

![Impact of electrical crosstalk in the detectors' readout electronics in CMB observations. For each Stokes parameter (top to bottom row) we show the input sky map (left), the sky map reconstructed from TOD affected by crosstalk (middle), and their difference (right).\label{fig:crosstalk}](crosstalk_zoom.png)

In the absence of instrumental systematic effects, the TOD of a single detector $d$ is modeled as:

\begin{equation}
d_t  = g_t[I^{\rm obs}(\hat{\mathbf{n}}_t)+\cos2\psi_t Q^{\rm obs}(\hat{\mathbf{n}}_t) + \sin2\psi_t U^{\rm obs}(\hat{\mathbf{n}}_t)] +n_t,
\label{eq:datamodel}
\end{equation}

where a $t$ subscript denotes a time-dependent quantity. $\mathbf{n}$ is the vector describing the detector pointing and $\psi$ is the polarization angle that describes the effective orientation of a polarization sensitive detector with respect to the input sky coordinate system. This depends on the orientation of the detector itself but also on the orientation of optical elements that modulate the incoming polarized signal (e.g. HWPs). In absence of systematics, the calibration factor (gain) is $g_t=1$ to preserve the calibration of the input signal. The exact value of the pixelized input Stokes parameters at $\mathbf{n}_t$ is determined through nearest grid point interpolation. White or correlated noise $n_t$ can be added to TOD according to the instrument specifications.

### Instrumental systematic effects

The instrumental systematic effects implemented in `s4cmb` are available in the `systematics.py` module. This includes electrical crosstalk in the detectors' readout electronics, gain misestimation or drifting of their values in time, pointing errors, distortions of the beam shape compared to the shape of the expected beam $B$ used to convolve the input sky and misestimation of their position in the focal plane, errors in the polarization angle estimation. These effectively modify on the fly the $\psi_t$, $\mathbf{n}_t$, $g_t$ in \autoref{eq:datamodel} compared to their expected value determined by the instrument design and scanning strategy when creating the TODs. Effects of beam distortions are conversely modelled through a Taylor expansion of $B$. Further details of the mathematical modelling of systematics are given in @mirmelstein2020 and @bicep-beam. We encourage users to add more effects and integrate their work in the package through pull requests on GitHub.

# Bootcamp
We release a [bootcamp](https://github.com/JulienPeloton/s4cmb-resources) dedicated to the package in two parts (beginners and advanced users) that include notebooks describing the basic parts of the API, and providing ready-to-use examples for the major code functionalities. A notebook to create \autoref{fig:crosstalk} can be found [here](https://github.com/JulienPeloton/s4cmb-resources/blob/master/Part1/s4cmb_crosstalk_05-joss.ipynb).

# Acknowledgements

We thank Neil Goeckner Wald for contributing to a large part of the scanning
strategy module, and the Polarbear collaboration for fostering the developments that led to this package. We acknowledge support from the European Research Council under the European Union’s Seventh Framework Programme (FP/2007-2013) / ERC Grant Agreement No. [616170]. GF acknowledges the support of the European Research Council under the Marie Skłodowska Curie actions through the Individual Global Fellowship No. 892401 PiCOGAMBAS.

# References
