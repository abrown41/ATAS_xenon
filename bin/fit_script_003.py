# -*- coding: utf-8 -*-
"""
Illustration of fit routine used for the evaluation of my measurement data.
@author: mhart

ACB: the csv files to be read by this script can be generated from the RMT
TA_spect_len_0001_z files using the makedf.py utility. The csv files should then
be named int1.3.csv (for intensity 1.3) and all stored in the same directory.
This script can then be run in that directory to do the fitting.
"""

#-----------------------------------------------------------------------
#                    Imports
#-----------------------------------------------------------------------
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as constants
import pandas as pd
import sys
import matplotlib.pyplot as plt



#-----------------------------------------------------------------------
#                   Pseudo Code
#-----------------------------------------------------------------------
def get_energy_axis():
    your_energy_axis = np.linspace(51.012496,67.000271,1962)
    return your_energy_axis

def getOD(intensity):
    """ Returns a whole time delay scan at a given intensity.
        Shape (td_axis_size, energy_axis_size). """
    fname = "OD.csv"
    df = pd.read_csv(fname)
    time_delay_axis = np.array([str(x) for x in df.columns[1:]])
    return df, time_delay_axis #df.transpose().values, time_delay_axis

def get_intensities():
    # return [1.3,1.6,1.9,2.2,2.5]
    return [2.2]

def get_roi(energy_axis,erange=(55.15,57.45)):
    for i,e in enumerate(energy_axis): 
        if e>erange[0]:
            first = i
            break
    for i,e in enumerate(energy_axis): 
        if e>erange[1]:
            last = i
            break

    return (first,last)
#-----------------------------------------------------------------------
#                   Constants
#-----------------------------------------------------------------------
path_length_density_product = 3e16 # in cm^-2  <-- estimated by OD and cross-section of T1 lines at iris 55
pldp_au = 0.77              # path-length-density-product in atomic units
alpha = constants['fine-structure constant'][0]
lineshape_constant = pldp_au/np.log(10)*4*np.pi*alpha
gamma_xe1 = 0.122 # from literature (Anderson 2001 I think?)

# resonance energies after calibration
# (retrieved by fitting a Lorentzian to the spectra far out of temporal overlap)
e_res = [55.38, 55.98, 57.27]

# set region of interest for fit
roi = slice(*get_roi(get_energy_axis()))
photonenergy = get_energy_axis()[roi]
#print (photonenergy)

# list of your intensity values
intensities = get_intensities()

#-----------------------------------------------------------------------
#                   Functions
#-----------------------------------------------------------------------
def wrap(phase, offset=0):
    """ Opposite of np.unwrap. Restrict phase to [-2*pi, 2*pi]. """
    return ( phase + np.pi + offset) % (2 * np.pi ) - np.pi - offset

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
    lineshape = (gamma/2*np.cos(phi) - (energy_axis-resonance_energy)*np.sin(phi)) / ((energy_axis-resonance_energy)**2 + gamma**2/4)
    return  z * lineshape

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
    # for e, energy in enumerate(e_res):
    for e, energy in [(2,e_res[2]), (1,e_res[1]), (0, e_res[0])]:
        model += DCM_lineshape(energy_axis, z[e]*gamma[e], phi[e], energy, gamma[e])
        model *= energy_axis* lineshape_constant
        model += params[-1] # add non-resonant background
    return model

#-----------------------------------------------------------------------
#                   Fit setup
#-----------------------------------------------------------------------

lower_bounds = [1e-6,  -2*np.pi, 0.9*gamma_xe1]*3 + [-15]
upper_bounds = [np.inf, 2*np.pi, 1.1*gamma_xe1]*3 + [20]
bounds = (lower_bounds, upper_bounds)

p_initial = [1, 0, gamma_xe1] + [0.1, 0, gamma_xe1]  + [0.01, 0, gamma_xe1]  + [0]
# initial_model = fit_lineshapes(photonenergy, *p_initial)
optimum = []
fit_errs = []
td = []

#-----------------------------------------------------------------------
#                   Fit loop
#-----------------------------------------------------------------------

df=pd.DataFrame()

intensity=intensities
OD,time_delay_axis = getOD(intensity)

for time_delay in time_delay_axis:
    freq_marginal=OD[str(time_delay)]
    plt.plot(photonenergy, freq_marginal[roi], label='Intensity Average')
    p_init, pcov = curve_fit(fit_lineshapes, photonenergy, freq_marginal[roi], p_initial, maxfev=500000, bounds=bounds)
    # plt.plot(photonenergy,fit_lineshapes(photonenergy,*p_initial),'r', label='fit init')
    plt.plot(photonenergy,fit_lineshapes(photonenergy,*p_init),'g', label='fit')
    
    plt.legend()
    plt.show()

z0=p_init[:-1:3]
print(z0)
