# -*- coding: utf-8 -*-
"""
Uses dipole.csv and field.csv (output from makedf.py) to calculate the OD
at each time delay in a time delay scan, perform the fitting procedure and
save the time delay scan, the fits at each time delay and the fitting parameters
to .csv files named OD.csv, OD_fit.csv and fit_params.csv respectively.

In calculating the OD, the dipole is tapered with an exponential decay to simulate
Auger Decay.

Adapted from a script by Max.
"""

#-----------------------------------------------------------------------
#                    Imports
#-----------------------------------------------------------------------
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as constants
import pandas as pd
import matplotlib.pyplot as plt
import helper_functions as hf


#-----------------------------------------------------------------------
#                   Constants
#-----------------------------------------------------------------------
gamma_xe1 = 0.122 # from literature (Anderson 2001 I think?)
alpha = constants['fine-structure constant'][0]

# set region of interest for fit
roi = slice(*hf.get_roi(hf.get_energy_axis(),erange=(60,62)))
photonenergy = hf.get_energy_axis()[roi]

#-----------------------------------------------------------------------
#                   Fit setup
#-----------------------------------------------------------------------
lower_bounds = [1e-6,  -np.pi, 0.5*gamma_xe1] + [-1]
upper_bounds = [np.inf, np.pi, 2*gamma_xe1] + [1]
bounds = (lower_bounds, upper_bounds)

p_init = np.array([1, 0, gamma_xe1] + [0])
params = []
errors = []
#-----------------------------------------------------------------------
#                   Read Data 
#-----------------------------------------------------------------------

df1 = pd.read_csv("dipole.csv")
df2 = pd.read_csv("field.csv")

#-----------------------------------------------------------------------
#                   Delay loop 
#-----------------------------------------------------------------------

# Get delays
delay_strings = df1.columns.values[1:]
td = delay_strings.astype(float)
OD_sim = []

# prepad signal with zeros to improve resolution
file_length=len(df1)
slen = 2**21
nzero = slen - file_length

for delay in delay_strings:
    t_zero = 1 + float(delay) + 5.75

    xuv_dipole = df1[delay]
    time = hf.au_to_fs(df1['Time']) - t_zero
    efield = df2[delay]

#-----------------------------------------------------------------------
#                   OD with exponential decay
#-----------------------------------------------------------------------

    # decay setup
    lifetime = hf.au_to_fs(1/(hf.ev_to_au(gamma_xe1)/2))
    decay = np.exp(- time/lifetime)
    decay[time<0] = 1

    # blackman window
    # window = np.blackman(file_length)

    # prepad to improve FT resolution
    prepad = np.zeros(nzero)
    field      = np.concatenate((prepad,  efield))
    dipole_pad = np.concatenate((prepad, decay * xuv_dipole))

    # FT to frequency domain
    dw = np.fft.fft(dipole_pad)
    ew = np.fft.fft(field)
    
    t_mult = time[1] - time[0]
    w = np.fft.fftfreq(field.shape[0], hf.fs_to_au(t_mult)/(2*np.pi))
    w_ev = hf.au_to_ev(w)
    w_roi = w_ev[np.all([(w_ev > 59), (w_ev < 62)], axis=0)]
    
    OD = -4*np.pi*alpha*w*np.imag(dw/ew)
    OD_roi = OD[np.all([(w_ev > 59), (w_ev < 62)], axis=0)]
    OD_sim.append(OD_roi)

#-----------------------------------------------------------------------
#                   Fit 
#-----------------------------------------------------------------------
    
    # restrict fit parameters to resonable values
    p_init[:-1:3] = abs(p_init[:-1:3])
    p_init[2:-1:3] = abs(p_init[2:-1:3])    

    popt, pcov = curve_fit(hf.fit_lineshapes, w_roi, OD_roi, 
                             p_init, maxfev=1000000, bounds=bounds)

    params.append(popt)
    errors.append(np.sqrt(np.abs(np.diag(pcov))))
    
    # reuse fit result in next iteration:
    # p_init = popt
    
params = np.array(params)
errors = np.array(errors)

#%%
#-----------------------------------------------------------------------
#                   Plots
#-----------------------------------------------------------------------
td*=-1
fig, (ax1,ax2,ax3) = plt.subplots(1, 3)

ax1.plot(td, params[:,0])
ax1.set_title('Line strength')

ax2.plot(td, params[:,1])
ax2.set_title('Phase (rad)')
ax2.set_xlabel('Time delay (fs)')

ax3.plot(td, params[:,2])
ax3.set_title('Line width (eV)')

plt.figure(3)
plt.plot(td, params[:,3])
plt.ylabel('Background')
plt.xlabel('Time delay (fs)')

#-----------------------------------------------------------------------
#                   Image
#-----------------------------------------------------------------------

paramdict = {'cmap':'jet', 'shading':'nearest'}
fontdict = {'color':'w', 'fontsize':10}

OD_sim = np.array(OD_sim)
OD_fit = np.zeros((params.shape[0], w_roi.size))
for p, popt in enumerate(params):
    OD_fit[p] = hf.fit_lineshapes(w_roi, *popt)


fig, ax = plt.subplots(nrows=1,ncols=2, num=5)
fig.subplots_adjust(right=0.9, left=0.1, top=0.9, bottom=0.15, wspace=0.2)

im = ax[0].pcolor(w_roi-5.5, td, OD_sim, **paramdict, vmin=-0.04, vmax=0.25)
vmin, vmax = im.get_clim()
im2 = ax[1].pcolor(w_roi-5.5, td, OD_fit, **paramdict, vmin=vmin, vmax=vmax)

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
cbar = fig.colorbar(im2, cax=cbar_ax)
cbar_ax.set_title('OD')
cbar.set_ticks([0,0.1,0.2])

ax[0].set_ylabel(r'Time delay (fs)')
ax[0].set_xlabel(r'Photon energy (eV)')
ax[1].set_xlabel(r'Photon energy (eV)')

ax[0].set_title('Simulation', fontsize=12, pad=1)    
ax[1].set_title('Fit', fontsize=12, pad=1)  

plt.show()

#-----------------------------------------------------------------------
#                   Save Time Delay Scan + Fit
#-----------------------------------------------------------------------

OD_df = pd.DataFrame()
OD_df['Energy'] = w_roi
for col in range(len(OD_sim[:,0])):
    delay = str(td[col])
    OD_df[delay] = OD_sim[col,:]

OD_df.to_csv('OD.csv', index=False)

fit_df = pd.DataFrame()
fit_df['Energy'] = w_roi
for col in range(len(OD_fit[:,0])):
    delay = str(td[col])
    fit_df[delay] = OD_fit[col,:]

fit_df.to_csv('OD_fit.csv', index=False)

#-----------------------------------------------------------------------
#                   Save Fitting Parameters
#-----------------------------------------------------------------------

params_df = pd.DataFrame()
params_df['Time Delays'] = td 
params_df['Line Strength'] = params[:,0]
params_df['Phase'] = params[:,1]
params_df['Line Width'] = params[:,2]

params_df.to_csv('fit_params.csv', index=False)
