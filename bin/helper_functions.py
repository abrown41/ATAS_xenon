# -*- coding: utf-8 -*-

"""

Functions to help with the analysis



Created on Thu Nov  4 10:00:07 2021



@author: mhart

"""

import numpy as np
import pandas as pd
import scipy.constants as cnt
import matplotlib.pyplot as plt
import scipy.optimize as sp
import glob
from argparse import ArgumentParser as AP
from pathlib import Path


# -----------------------------------------------------------------------
#                   Constants
# -----------------------------------------------------------------------

pldp_au = 0.77              # path-length-density-product in atomic units
alpha = cnt.constants.fine_structure
lineshape_constant = pldp_au/np.log(10)*4*np.pi*alpha

# -----------------------------------------------------------------------
#                   Functions
# -----------------------------------------------------------------------


def au_to_fs(time_in_au):
    """
    converts atomic unit to femto seconds

    :param time_in_au: time to convert

    :type time_in_au: float

    :return: time in femto seconds

    :rtype: float

    """
    return time_in_au * cnt.value("atomic unit of time") / cnt.femto


def fs_to_au(time_in_femto_seconds):
    """
    converts femto seconds to atomic unit

    :param time_in_femto_seconds: time to convert

    :type time_in_femto_seconds: float

    :return: time in atomic units

    :rtype: float

    """
    return time_in_femto_seconds * cnt.femto / cnt.value("atomic unit of time")


def au_to_ev(energy_in_au):
    """
    converts atomic units to electron volt



    :param energy_in_au: energy to convert

    :type energy_in_au: float

    :return: returns energy in electron volt

    :rtype: float

    """
    ev_to_joule = energy_in_au * cnt.value("atomic unit of energy")
    return ev_to_joule / cnt.value("electron volt-joule relationship")


def ev_to_au(energy_in_ev):
    """
    converts electron volts to atomic units

    :param energy_in_ev: energy to convert

    :type energy_in_ev: float

    :return: returns energy in atomic units

    :rtype: float
    """
    ev_to_joule = energy_in_ev * cnt.value("electron volt-joule relationship")
    return ev_to_joule / cnt.value("atomic unit of energy")


def cos2_window(index1, index2, width, size):

    indmax = np.maximum(index1, index2)
    indmin = np.minimum(index1, index2)

    arr1 = np.sin(np.pi/2/width*np.array(range(width)))**2.
    arr2 = np.ones((indmax-indmin))
    arr3 = np.cos(np.pi/2/width*np.array(range(width)))**2.

    checksum = size-indmax+indmin-2*width
    arr4 = np.zeros((checksum))

    if np.greater_equal(checksum, 0.):
        return np.roll(np.concatenate((arr1, arr2, arr3, arr4)), indmin-width)

    else:
        print('Cos2 Window error. Sizes not okay.')
        return np.zeros((size))


def taper_cos2(width, array):

    number_of_zeros = int(width / 100.)
    width -= number_of_zeros
    modulation = np.ones(array.shape[0])

    if not width == 0:
        taper_mod = np.arange(1, width + 1)
        taper_mod = np.cos(np.pi / (2.0 * width) * taper_mod) ** 2
        modulation[array.shape[0] - (width + number_of_zeros)
                                     :array.shape[0] - number_of_zeros] = taper_mod
        modulation[number_of_zeros:width + number_of_zeros] = taper_mod[::-1]

    modulation[array.shape[0] - number_of_zeros:] = 0.0
    modulation[0:number_of_zeros] = 0.0
    return modulation


def TDDM_Reconstruction(data, ax_off, ax_mult):

    if len(data.shape) == 1:

        p_data = np.concatenate((np.zeros(int(ax_off/ax_mult)), data))

        axis = au_to_fs(np.fft.fftfreq(
            p_data.shape[0], ev_to_au(ax_mult)/(2*np.pi)))

        tddm = np.fft.fft(p_data, norm='ortho')

        return (tddm, axis)

    if len(data.shape) == 2:

        zeros_front = np.zeros((int(ax_off/ax_mult), data.shape[1]))

        p_data = np.concatenate((zeros_front, data), axis=0)

        axis = au_to_fs(np.fft.rfftfreq(
            p_data.shape[0], ev_to_au(ax_mult)/(2*np.pi)))

        tddm = np.zeros_like(p_data, dtype=complex)

        tddm = np.fft.fft(p_data, axis=0, norm='ortho')

        return (tddm, axis)


def get_energy_axis():

    your_energy_axis = np.linspace(51.012496, 67.000271, 1962)

    return your_energy_axis


def getOD(intensity):
    """ Returns a whole time delay scan at a given intensity.

        Shape (td_axis_size, energy_axis_size). """

    fname = "int"+str(intensity)+".csv"

    df = pd.read_csv(fname)

    time_delay_axis = np.array([float(x) for x in df.columns])

    return df.transpose().values, time_delay_axis


def get_intensities():

    return [1.3, 1.6, 1.9, 2.2, 2.5]


def get_roi(energy_axis, erange=(60, 62)):

    for i, e in enumerate(energy_axis):

        if e > erange[0]:

            first = i

            break

    for i, e in enumerate(energy_axis):

        if e > erange[1]:

            last = i

            break

    return (first, last)


def wrap(phase, offset=0):
    """ Opposite of np.unwrap. Restrict phase to [-2*pi, 2*pi]. """

    return (phase + np.pi + offset) % (2 * np.pi) - np.pi - offset


def DCM_lineshape(energy_axis, z, phi, resonance_energy, gamma):
    """

    Dipole control model (DCM) line shape function for a single absorption line



    Parameters

    ----------

    energy_axis : the array of values that defines the photon energy axis

    z : line strength

    phi : dipole phase

    resonance_energy : resonance energy of the absorption line

    gamma : line width



    Returns

    -------

    np.array size of energy axis

        line shape function as a function of photon energy

    """

    lineshape = (gamma/2*np.cos(phi) - (energy_axis-resonance_energy)
                 * np.sin(phi)) / ((energy_axis-resonance_energy)**2 + gamma**2/4)

    return z * lineshape


# resonance energies after calibration

# (retrieved by fitting a Lorentzian to the spectra far out of temporal overlap)
e_res = [60.88]


def fit_lineshapes(energy_axis, *params):
    """

    Fit function to extract line shape parameters from several absorption lines

    from the measurement data. I omitted the convolution with the experimental

    spectral resolution which only marginally affects the fit results anyway.



    Parameters

    ----------

    energy_axis : the array of values that defines the photon energy axis

    *params : list of fit parameters, size: 3*N + 1 where N is the number of lines



    Returns

    -------

    model : np.array size of energy axis

        Calculates an optical density as a function of photon energy by adding

        up the line shape functions of N absorption lines. Includes a constant

        offset to fit the non-resonant background.

    """

    model = np.zeros(energy_axis.shape)

    z = params[:-1:3]

    phi = params[1:-1:3]

    gamma = params[2:-1:3]

    for e, energy in enumerate(e_res):

        model += DCM_lineshape(energy_axis,
                               z[e]*gamma[e], phi[e], energy, gamma[e])

        model *= energy_axis * lineshape_constant

        model += params[-1]

    return model


def read_command_line():
    parser = AP()
    parser.add_argument('-p', '--plot', help="show a plot",
                        action='store_true', default=False)
    parser.add_argument('-i', '--IR_intensity',
                        type=float, help="IR intensity")
    return vars(parser.parse_args())


def moving_average(input_data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(input_data, window, 'same')


def smooth_data(wdata, nits=8):
    out_data = wdata
    for ii in np.arange(nits):
        out_data = moving_average(out_data, 4)
    return out_data


def strip_num(filename):
    """
    Obtain the time delay from the directory name
    """
    import re
    pat = re.compile("[0-9]{1}.[0-9]{2}")
    m = pat.search(filename)
    pat = re.compile("[+,-]")
    n = pat.search(filename)
    return(float(n.group(0)+m.group(0)))


def order_files(filelist):
    """
    Sort the files into increasing order of time delay
    """
    from operator import itemgetter
    td_list = []
    for fname in filelist:
        td_list.append((fname, strip_num(fname)))

    return(sorted(td_list, key=itemgetter(1), reverse=True))


def grab_data(fname):
    df = pd.read_csv(fname, delim_whitespace=True)
    return (df)
