"""
This script to removes IR contribution from the IR+XUV dipole to give only the XUV contribution.
Uses IR_only_dipole.csv and .csv file containing IR+XUV dipole at each time delay (generated using makedf.py).
IR+XUV dipole .csv file must be named after IR intensity (e.g. int1.3.csv)
Run with:
    python remove_IR.py <path_to_IR_file> <path_to_IR+XUV_file>
I *think* the output can then be used with makeOD.py to get TAS. Need to check...
"""

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser as AP
import pandas as pd
from pathlib import Path

def read_command_line():
    parser = AP()
    parser.add_argument('IR_only_file',help="path to .csv file containing IR_only dipole at each IR intensity")
    parser.add_argument('XUV_file',help="path to .csv file containing IR+XUV dipole at each time delay, named after IR intensity e.g. int1.3.csv (generated using makedf.py)")
    return vars(parser.parse_args())

def strip_num(filename):
    "Gets IR intensity from filename"
    import re
    pat=re.compile("[0-9]{1}.[0-9]{1}")
    m=pat.search(filename)
    return(float(m.group(0)))

args=read_command_line()

IR_df = pd.read_csv(Path(args["IR_only_file"]))
IR_intensity = str(strip_num(args["XUV_file"]))
IR_dipole = IR_df[IR_intensity]

df = pd.read_csv(Path(args["XUV_file"]))

# if len(IR_df) != len(df):
#     print('Attempting to subtract two files of different lengths. Check Final time of calculation is consistent.')

ndf = pd.DataFrame()
ndf[0] = df["Time"]
for col in range(len(df.columns)-1):
    ndf[col+1] = df.iloc[:,col+1] - IR_dipole

ndf.to_csv("XUV_only_"+IR_intensity+".csv", index=False, sep=',', header=list(df))