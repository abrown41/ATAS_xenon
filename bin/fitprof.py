"""
Original Script for fitting

This script extracts the time-dependent oscillations in the absorption peaks and
fits them with a fano/lorentzian profile so that the phase, amplitude and
linewidth can be extracted. Run with python fitprof.py <list of rmt directories>
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

    return(sorted(td_list, key=itemgetter(1)))

def grabData(fname,erange=(55.15,55.62)):
    df = pd.read_csv(fname, delim_whitespace=True, header=None)
    first = (min(df[df[0]>erange[0]].index.tolist()))
    last = (min(df[df[0]>erange[1]].index.tolist()))
    df=df.truncate(before=first,after=last)
    df[1] *= -1
    return (df)

    

def fitfunc(w,t,wT,gam,z,phi,C):
    tline = 0.5*gam *np.cos(phi) + (w-wT) * np.sin(phi)
    bline = (w-wT)**2 + gam*gam/4
    prefac = 0.77 * 4*np.pi * w/ (137 * np.log(10))
    return C + prefac *  z  * tline / bline

def getFit(df,t,wT):
    try :
#                [  t ,       wT,   ,  gam,z, phi, C]        [  t ,       wT,   ,  gam,z,phi,C] 
        bounds = ([t-0.00001,wT-0.001,0.090,0, -np.pi, -1 ], [ t+0.00001, wT+0.001,0.110,1, 1, 1 ])
        popt,pcov = sp.curve_fit(fitfunc,df[0],df[1],bounds=bounds,method='trf')
        return(popt)
    except:
        return np.array([-999, -999, -999, -999, -999, -999])

#wT=55.376# energy of resonance (eV)
#erange=(55.15,55.62)
#wT=55.98 # energy of resonance (eV)
#erange=(55.84,56.16)
wT=57.27
erange=(57.10,57.45)
fs=41.34 # conversion factor from a.u to fs.

args = read_command_line()
filelist=args["files"]
filelist= order_files(filelist)

z0=[]
phi=[]
gam=[]
times = []
fig, (ax1,ax2,ax3) = plt.subplots(1, 3)

for f,t in filelist:
    df = grabData(f+"/TA_spect_len_0001_z",erange=erange)
    popt = getFit(df,t*fs,wT)
    if not np.any(popt < -900) :
        z0.append(popt[3])
        phi.append(popt[4])
        gam.append(popt[2])
        times.append(t)

ndf=pd.DataFrame()
ndf["t"] = times
ndf["z0"] = z0
ndf["phi"] = phi
ndf["gam"] = gam

ndf.to_csv(args["output"])
if args["plot"]:
    ax1.plot(times,z0)
    ax2.plot(times,phi)
    ax3.plot(times,gam)
    plt.show()

#fit = fitfunc(df[0], popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
#plt.plot(df[0],df[1],label="data")
#plt.plot(df[0],fit,label="fit")
#plt.legend()
#plt.show()

