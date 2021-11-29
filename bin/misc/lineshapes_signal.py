"""
LH: this script reads in the OD and
    transforms the line profile back into a signal.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as constants
import pandas as pd
import sys
import matplotlib.pyplot as plt


def grabData(fname,erange=(55.15,55.62)):
    df = pd.read_csv(fname)#, delim_whitespace=True, header=None)
    first = (min(df[df['Energy']>erange[0]].index.tolist()))
    last = (min(df[df['Energy']>erange[1]].index.tolist()))
    df=df.truncate(before=first,after=last)
    # df['05.00'] *= -1
    return (df)


wT=55.376# energy of resonance (eV)
erange=(54.879,55.71)
df = grabData('OD.csv', erange=erange)

slen = 2**21
nzero = slen - len(df['05.00'])
# convolve the signal with a blackman window to smooth it
# window = np.blackman(len(dipole))
prepad = np.zeros(nzero)



plt.plot(df['Energy'], df['05.00'])
plt.show()
line = np.concatenate((prepad, df['05.00']))
signal = np.fft.fft(line)

plt.plot(signal)
plt.show()
