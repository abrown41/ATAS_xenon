"""
Read dipole and field data from RMT calculations, and compute the OD
Output into csv file with each column the OD for a given time delay
"""

import pandas as pd
import numpy as np
from scipy.constants import physical_constants as constants

alpha = constants['fine-structure constant'][0]

# load the dipole and field data from file
dipole = pd.read_csv("dipole.csv")
field = pd.read_csv("field.csv")


# prepad signal with zeros to improve resolution
slen = 2**21
nzero = slen - len(dipole)
prepad = np.zeros(nzero)

# convolve the signal with a blackman window to smooth it
window = np.blackman(slen)

# calculate the energy axis (in a.u.)
delta_t = 0.01  # time step in a.u.
freqs=np.arange((2**18))*(((2.0*np.pi)/(2**19))/delta_t)

# set the region of interest
roi = slice(6257,8219)
freqs=freqs[roi]

# dataframe to hold the output data 
df = pd.DataFrame() 

for delay in dipole:
# prepad & window
    dt = window * np.concatenate((prepad, dipole[delay]))
    et = window * np.concatenate((prepad, field[delay]))
# FFT and select slice
    dw = np.fft.fft(dt)[roi]
    ew = np.fft.fft(et)[roi]
# Calculate OD
    OD = -4*np.pi*alpha*freqs*np.imag(dw/ew)
    df[delay] = OD

df.to_csv("OD.csv",index=False)
