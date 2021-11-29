"""
LRH: This script reads in the IR-only dipole and the IR+XUV dipole, removes the IR contribution
    from the dipole and generates the optical density with only the XUV contribution (OD.csv) which
    can then be used to fit the lineshapes. Run from directory containing the dipole files.
    Edits to lines 15-20 will change the files read in.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import physical_constants as constants

alpha = constants['fine-structure constant'][0]

df = pd.read_csv("IR_only_dipole.csv")
dIR = df["2.2"].values

field = pd.read_csv("EField.Xe_ATAS00002200", delim_whitespace=True)
time = field['Time'].values
dtot = pd.read_csv('coherent.txt', delim_whitespace=True)['0001_z'].values

dipole = dtot-dIR

cutoff=150000
dipole = dipole[:cutoff]
field=field[:cutoff]

# prepad signal with zeros to improve resolution
slen = 2**21
nzero = slen - len(dipole)
# convolve the signal with a blackman window to smooth it
window = np.blackman(len(dipole))
prepad = np.zeros(nzero)


# calculate the energy axis (in a.u.)
delta_t = 0.01  # time step in a.u.
freqs=np.arange((2**20))*(((2.0*np.pi)/(2**21))/delta_t)

# set the region of interest
roi = slice(6257,8219)
freqs=freqs[roi]
print (freqs[0]*27.212, freqs[-1]*27.212)

# dataframe to hold the output data 
df = pd.DataFrame() 

dt =  np.concatenate((prepad, window * dipole))
et =  np.concatenate((prepad,  field['0001_z'].values))

# FFT and select slice
dw = np.fft.fft(dt)[roi]
ew = np.fft.fft(et)[roi]

# Calculate OD
OD = -4*np.pi*alpha*freqs*np.imag(dw/ew)

df['Energy'] = freqs*27.212
df["05.00"] = OD

plt.plot(freqs*27.212,OD)
plt.show()
df.to_csv("OD.csv",index=False)


