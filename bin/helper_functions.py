# -*- coding: utf-8 -*-

"""

Functions to help with the analysis



Created on Thu Nov  4 10:00:07 2021



@author: mhart

"""

import numpy as np
import scipy.constants as cnt
from argparse import ArgumentParser as AP

# -----------------------------------------------------------------------
#                   Constants
# -----------------------------------------------------------------------

pldp_au = 0.77              # path-length-density-product in atomic units
alpha = cnt.constants.fine_structure
lineshape_constant = pldp_au/np.log(10)*4*np.pi*alpha

# resonance energies after calibration

# (retrieved by fitting a Lorentzian to the spectra far out of temporal overlap)

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
    global e_res

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
    global e_res
    parser = AP()
    parser.add_argument('-p', '--plot', help="show the plotted data",
                        action='store_true', default=False)
    parser.add_argument('-o', '--output', help="save data to file",
                        action='store_true', default=True)
    parser.add_argument('-i', '--IR_intensity',
                        type=float, help="IR intensity")
    parser.add_argument('--energy_shift',
                        help="energy in eV by which to shift the energy axis",
                        type=float, default=-5.5)
    parser.add_argument('--e_res',
                        help="positions of resonances in eV (before shift)",
                        type=list, default=[60.88])

    args = vars(parser.parse_args())

    e_res = args['e_res']
    roi_lo = min(e_res) - 1.5
    roi_hi = max(e_res) + 1.5
    args['roi'] = [roi_lo, roi_hi]
    e_res = [e + args['energy_shift'] for e in e_res]
    return args


def moving_average(input_data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(input_data, window, 'same')


def smooth_data(wdata, nits=8):
    out_data = wdata
    for ii in np.arange(nits):
        out_data = moving_average(out_data, 4)
    return out_data
