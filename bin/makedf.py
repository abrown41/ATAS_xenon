"""
This script reads the TA_spect_len_0001_z files from the individual directories,
and puts them together into a dataframe, which is output to a csv file. This csv
file can then be read by the fit_script utility that max sent. 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as sp
from argparse import ArgumentParser as AP

def read_command_line():
    parser = AP()
    parser.add_argument('files', nargs="+",help="list of directories containing output for plotting")
    parser.add_argument('-o','--output', type=str,help="output file for fit parameters, default is fitparams.csv",default="fitparams.csv")
    parser.add_argument('-p','--plot', help="show a plot of the three parameters",action='store_true',default=False)
    return vars(parser.parse_args())


def strip_num(filename):
    import re
    pat=re.compile("[0-9]{2}.[0-9]{3}")
    m=pat.search(filename)
    pat=re.compile("[+,-]")
    n=pat.search(filename)
    return(float(n.group(0)+m.group(0)))


def order_files(filelist):
    from operator import itemgetter
    td_list=[]
    for fname in filelist:
        td_list.append((fname,strip_num(fname))) 

    return(sorted(td_list, key=itemgetter(1), reverse=True))

def grabData(fname):
    df = pd.read_csv(fname, delim_whitespace=True, header=None)
    df[1] *= -1
    return (df)

args = read_command_line()
filelist=args["files"]
filelist= order_files(filelist)

ndf = pd.DataFrame()
for f,t in filelist:
    df = grabData(f+"/TA_spect_len_0001_z")
    ndf[str(t)] = df[1]

#ndf.insert(0, "Energy", df[0])
ndf.to_csv(args["output"],index=False)


