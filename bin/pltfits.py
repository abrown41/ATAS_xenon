"""
Script to plot specific parameters from the output of fitprof.py"
"""

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser as AP
import pandas as pd
from pathlib import Path

def read_command_line():
    parser = AP()
    parser.add_argument('dirs', nargs="+",help="list of directories containing output for plotting")
    parser.add_argument('-f','--file',help="filename of fitting parameters csv, default= 'fitparmams.csv'", default='fitparams.csv')
    parser.add_argument('-e','--exp', nargs="+",help="csv files containing experimental data",default=None)
    parser.add_argument('-p','--param', type=str,help="parameter for plotting: default is z0 (options phi, gam)",default="z0")
    return vars(parser.parse_args())

args=read_command_line()

for f in args["dirs"]:
    p=Path(f)
    df = pd.read_csv(p / args["file"])
    plt.plot(df["t"],df[args["param"]],label=f)

#count = 0 
#if args["exp"]:
#    for f in args["exp"]:
#        df = pd.read_csv(f,header=None)
#        plt.plot(df[0],df[1],cols[count]+".",label="exp"+f.split("/")[-1][:-4])
#        count+=1

plt.legend()
plt.show()

