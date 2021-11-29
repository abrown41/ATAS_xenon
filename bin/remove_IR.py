"""
This script to removes IR contribution from the IR+XUV dipole to give only the XUV contribution.
Uses IR_only_dipole.csv and .csv file containing IR+XUV dipole at each time delay (generated using makedf.py).
IR+XUV dipole .csv file must be named after IR intensity (e.g. int1.3.csv)
Run with:
    python remove_IR.py <path_to_IR_file> <path_to_IR+XUV_file> <path_to_EField_file>
Creates two .csv files, one with XUV only dipole, and one with XUV only OD.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from argparse import ArgumentParser as AP
from scipy.constants import physical_constants as constants

def read_command_line():
    parser = AP()
    parser.add_argument('IR_only_file',help="path to .csv file containing IR_only dipole at each IR intensity")
    parser.add_argument('XUV_file',help="path to .csv file containing IR+XUV dipole at each time delay, named after IR intensity e.g. int1.3.csv")
    parser.add_argument('Efield_file',help="path to .csv file containing EField at each time delay, named after IR intensity e.g. EField1.3.csv")
    return vars(parser.parse_args())

def strip_num(filename):
    "Gets IR intensity from filename"
    import re
    pat=re.compile("[0-9]{1}.[0-9]{1}")
    m=pat.search(filename)
    return(float(m.group(0)))

args=read_command_line()

alpha = constants['fine-structure constant'][0]

# Get dipole/Efield files:
IR_df = pd.read_csv(Path(args["IR_only_file"]))
# IR_intensity = str(strip_num(args["XUV_file"]))
IR_intensity=str(2.2)
IR_dipole = IR_df[IR_intensity]

EField = pd.read_csv(Path(args["Efield_file"]))
df = pd.read_csv(Path(args["XUV_file"]))

# get XUV only dipole

XUV_only_dipole = pd.DataFrame()
XUV_only_dipole[0] = df["Time"]
for col in range(len(df.columns)-1):
    XUV_only_dipole[col+1] = df.iloc[:,col+1] - IR_dipole

XUV_only_dipole.to_csv("XUV_only_dipole"+IR_intensity+".csv", index=False, sep=',', header=list(df))

# Calculated XUV only OD:

# Truncate data after 150000 time steps:
cutoff=150000
XUV_only_dipole.truncate(before=1, after=cutoff)
EField.truncate(before=1, after=cutoff)

# prepad signal with zeros to improve resolution
slen = 2**21
nzero = slen - len(XUV_only_dipole)
# convolve the signal with a blackman window to smooth it
window = np.blackman(len(XUV_only_dipole))
prepad = np.zeros(nzero)

# calculate the energy axis (in a.u.)
delta_t = 0.01  # time step in a.u.
freqs=np.arange((2**20))*(((2.0*np.pi)/(2**21))/delta_t)

# set the region of interest
roi = slice(6257,8219)
freqs=freqs[roi]

OD_df = pd.DataFrame()
# OD_df['Energy'] = freqs*27.212


for col in range(len(XUV_only_dipole.columns)-1):
    # Prepad and window
    dt = np.concatenate((prepad, window * XUV_only_dipole.iloc[:,col+1]))
    et =  np.concatenate((prepad,  EField.iloc[:,col+1]))
    # FFT and select slice
    dw = np.fft.fft(dt)[roi]
    ew = np.fft.fft(et)[roi]
    # Calculate OD
    OD = -4*np.pi*alpha*freqs*np.imag(dw/ew)
    OD_df[col] = OD

# headers=list(df)
# headers[0] = 'Energy'

OD_df.to_csv("XUV_only_OD"+IR_intensity+".csv", index=False, sep=',', header=list(df)[1:])



