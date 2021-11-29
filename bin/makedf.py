"""
This script reads the TA_spect_len_0001_z files from the individual directories,
and puts them together into a dataframe, which is output to a csv file. This csv
file can then be read by the fit_script utility that max sent. 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as sp
import glob
from argparse import ArgumentParser as AP

def read_command_line():
    parser = AP()
    parser.add_argument('files', nargs="+",help="list of directories containing output for plotting")

    parser.add_argument('-o','--output', type=str,help="output file for fit parameters, default is output.csv",default="output.csv")
    return vars(parser.parse_args())


# def strip_num(filename):
#     import re
#     pat=re.compile("[0-9]{2}.[0-9]{3}")
#     m=pat.search(filename)
#     pat=re.compile("[+,-]")
#     n=pat.search(filename)
#     return(float(n.group(0)+m.group(0)))

# def strip_num_IR(filename):
#     import re
#     pat=re.compile("[0-9]{1}.[0-9]{1}")
#     m=pat.search(filename)
#     return(float(m.group(0)))

def strip_num_OR_dir(filename):
    import re
    pat=re.compile("[0-9]{2}.[0-9]{1}")
    m=pat.search(filename)
    return(float(m.group(0)))

def order_files(filelist):
    from operator import itemgetter
    td_list=[]
    for fname in filelist:
        # td_list.append((fname,'+05.00')) 
        td_list.append((fname,strip_num_OR_dir(fname))) 

    return(sorted(td_list, key=itemgetter(1), reverse=True))

def grabData(fname):
    df = pd.read_csv(fname, header=None)#, delim_whitespace=True, header=None)
    print(df)
    return (df)

args = read_command_line()
filelist=args["files"]
filelist= order_files(filelist)

ndf = pd.DataFrame()
for f,t in filelist:
    print(f)
    fname=glob.glob(f+"/OD.csv")[0]
    df = grabData(fname)
    ndf[str(t)] = df[1]

ndf.insert(0, "Energy", df[0])
ndf.drop(index=1)
ndf.to_csv(args["output"]+'.csv',index=False)


