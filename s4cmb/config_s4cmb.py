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

def boolise_it(entry):
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
    >>> config_dict = {'toto': 'True', 'tata': 'False', 'tutu': '1.2'}
    >>> boolise_it('True')
    True

    >>> boolise_it('False')
    False

    >>> boolise_it('1.2') #doctest: +NORMALIZE_WHITESPACE
    You assign 1.2 to bool, but is neither True or False...
    Return False by default, but I suggest that you you look at your ini file.
    False
    """
    if entry == 'False':
        return False
    elif entry == 'True':
        return True
    else:
        print('You assign {} to bool, but is neither True or '.format(entry) +
              'False... Return False by default, but I suggest that you ' +
              'you look at your ini file.')
        return False
    return out

def find_value_from_last_letter(string):
    """
    Routine to find the type of a variable in the parser.
    The variable should be entered in the format
    <name> = <value> <letter>
    where letter is used to determine the type of the variable:
    S(tring), F(loat), I(nteger), B(ool), N(one).
    """
    type_of_my_string = string[-1]
    if type_of_my_string == 'S':
        return string[:-2]
    elif type_of_my_string == 'F':
        return floatise_it(string[:-2])
    elif type_of_my_string == 'I':
        return intise_it(string[:-2])
    elif type_of_my_string == 'B':
        return boolise_it(string[:-2])
    elif type_of_my_string == 'N':
        return None

class NormaliseParser():
    """
    Class to handle s4cmb parsers.
    It converts initial dictionary into an object.
    """
    def __init__(self, config_dict):
        """
        Careful: keys should be in small caps (as returned by ConfigParser) to
        avoid confusion.

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
        >>> params = NormaliseParser(Config._sections['s4cmb'])
        >>> assert hasattr(params, 'nces')
        """
        for key, string_value in config_dict.iteritems():
            value = find_value_from_last_letter(string_value)
            setattr(self, key, value)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
