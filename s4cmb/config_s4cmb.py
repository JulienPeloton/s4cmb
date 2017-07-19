#!/usr/bin/python
"""
Module to normalise ini files containing parameter
values for running s4cmb in software mode.
Not required for the API.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""

def floatise_it(entry):
    """
    Convert a string to a float

    Parameters
    ----------
    entry : string
        String to convert

    Returns
    ----------
    float : float
        Converted entry.

    Examples
    ----------
    >>> print(floatise_it('1.2'))
    1.2
    """
    return float(entry)

def intise_it(entry):
    """
    Convert a string to a int

    Parameters
    ----------
    entry : string
        String to convert

    Returns
    ----------
    int : int
        Converted entry.

    Examples
    ----------
    >>> print(intise_it('10'))
    10
    """
    return int(entry)

def boolise_it(dic, entry):
    """
    Convert a string to a bool

    Parameters
    ----------
    entry : string
        String to convert. Has to be 'True' or 'False'.

    Returns
    ----------
    bool : bool
        Converted entry.

    Examples
    ----------
    >>> import numpy as np
    >>> config_dict = {'toto': 'True', 'tata': 'False'}
    >>> np.all(boolise_it(config_dict, 'toto'))
    True

    >>> np.all(boolise_it(config_dict, 'tata'))
    False

    >>> np.all(boolise_it(config_dict, 'titi'))
    False
    """
    if entry in dic:
        out = 'True' in dic[entry]
    else:
        out = False
    return out

class NormaliseS4cmbParser():
    """
    Class to handle s4cmb parser.
    It converts initial dictionary into an object.
    """
    def __init__(self, config_dict):
        """
        Careful: keys are in small caps (as returned by ConfigParser).

        Parameters
        ----------
        config_dict : dictionary
            dictionary coming from the ini file.
            Contains {key: val}, where val are strings
            (default to ConfigParser).

        Examples
        ----------
        >>> import ConfigParser
        >>> Config = ConfigParser.ConfigParser()
        >>> fns = Config.read('examples/simple_parameters.ini')
        >>> params = NormaliseS4cmbParser(Config._sections['s4cmb'])
        >>> assert hasattr(params, 'nCES')
        """

        ## Names
        self.input_filename = config_dict['input_filename']
        self.name_strategy = config_dict['name_strategy']
        self.pm_name = config_dict['pm_name']
        self.type_HWP = config_dict['type_hwp']
        self.ut1utc_fn = config_dict['ut1utc_fn']
        self.language = config_dict['language']
        self.start_date = config_dict['start_date']
        self.telescope_longitude = config_dict['telescope_longitude']
        self.telescope_latitude = config_dict['telescope_latitude']
        self.name_instrument = config_dict['name_instrument']

        ## Booleans
        self.do_pol = boolise_it(config_dict, 'do_pol')
        self.no_ileak = boolise_it(config_dict, 'no_ileak')
        self.no_quleak = boolise_it(config_dict, 'no_quleak')
        self.ext_map_gal = boolise_it(config_dict, 'ext_map_gal')
        self.verbose = boolise_it(config_dict, 'verbose')

        ## Integers
        self.ncrate = intise_it(config_dict['ncrate'])
        self.ndfmux_per_crate = intise_it(config_dict['ndfmux_per_crate'])
        self.nsquid_per_mux = intise_it(config_dict['nsquid_per_mux'])
        self.npair_per_squid = intise_it(config_dict['npair_per_squid'])
        self.beam_seed = intise_it(config_dict['beam_seed'])
        self.nCES = intise_it(config_dict['nces'])
        self.nside_out = intise_it(config_dict['nside_out'])

        ## Float
        self.FWHM_in = floatise_it(config_dict['fwhm_in'])
        self.FWHM = floatise_it(config_dict['fwhm'])
        self.fp_size = floatise_it(config_dict['fp_size'])
        self.projected_fp_size = floatise_it(
            config_dict['projected_fp_size'])
        self.freq_HWP = floatise_it(config_dict['freq_hwp'])
        self.angle_HWP = floatise_it(config_dict['angle_hwp'])
        self.telescope_elevation = floatise_it(
            config_dict['telescope_elevation'])
        self.sampling_freq = floatise_it(config_dict['sampling_freq'])
        self.sky_speed = floatise_it(config_dict['sky_speed'])
        self.width = floatise_it(config_dict['width'])

        ## Default None
        if config_dict['nside_in'] == 'None':
            self.nside_in = None
        else:
            self.nside_in = intise_it(config_dict['nside_in'])

        if config_dict['map_seed'] == 'None':
            self.map_seed = None
        else:
            self.map_seed = intise_it(config_dict['map_seed'])

class NormaliseXpureParser():
    """
    Class to handle xpure parser.
    It converts initial dictionary into an object.
    """
    def __init__(self, config_dict):
        """
        Careful: keys are in small caps (as returned by ConfigParser).

        Parameters
        ----------
        config_dict : dictionary
            dictionary coming from the ini file.
            Contains {key: val}, where val are strings
            (default to ConfigParser).

        Examples
        ----------
        >>> import ConfigParser
        >>> Config = ConfigParser.ConfigParser()
        >>> fns = Config.read('examples/xpure_parameters.ini')
        >>> params = NormaliseXpureParser(Config._sections['xpure'])
        >>> assert hasattr(params, 'node')
        """

        ## Names
        self.time = config_dict['time']
        self.queue = config_dict['queue']
        self.beam_file = config_dict['beam_file']
        self.bin_file = config_dict['bin_file']

        ## Integers
        self.node = intise_it(config_dict['node'])
        self.nproc_apo = intise_it(config_dict['nproc_apo'])
        self.nproc_scalar_to_spin = intise_it(
            config_dict['nproc_scalar_to_spin'])
        self.nproc_mll = intise_it(config_dict['nproc_mll'])
        self.nproc_xpure = intise_it(config_dict['nproc_xpure'])
        self.radius_apodization = intise_it(
            config_dict['radius_apodization'])
        self.lmax_user = intise_it(config_dict['lmax_user'])
        self.xpure_mode = intise_it(config_dict['xpure_mode'])
        self.fast = intise_it(config_dict['fast'])


if __name__ == "__main__":
    import doctest
    doctest.testmod()
