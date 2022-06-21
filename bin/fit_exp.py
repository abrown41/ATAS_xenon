# -*- coding: utf-8 -*-
"""
Script for fitting three peaks in experimental data.
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

args = hf.read_command_line()
intensity = args["IR_intensity"]

#-----------------------------------------------------------------------
#                   Constants
#-----------------------------------------------------------------------
gamma_xe1 = 0.122 # from literature (Anderson 2001 I think?)
alpha = constants['fine-structure constant'][0]
e_res = [55.38, 55.98, 57.30]

#-----------------------------------------------------------------------
#                   Fit setup
#-----------------------------------------------------------------------
lower_bounds = [1e-6,  -2*np.pi, 0.5*gamma_xe1]*3 + [-1]
upper_bounds = [np.inf, 2*np.pi, 2*gamma_xe1]*3 + [1]
bounds = (lower_bounds, upper_bounds)

p_init = np.array([1, 0, gamma_xe1] + [0.1, 0, gamma_xe1]  + [0.01, 0, gamma_xe1]  + [0])
params = []
errors = []
#-----------------------------------------------------------------------
#                   Data 
#-----------------------------------------------------------------------

df1 = pd.read_csv(f"OD{intensity}.csv")
w_ev = df1['Energy']

# set region of interest for fit
roi = [50,60]
w_roi = w_ev[np.all([(w_ev > roi[0]), (w_ev < roi[1])], axis=0)]

photonenergy = w_ev[roi]

#-----------------------------------------------------------------------
#                   Delay loop 
#-----------------------------------------------------------------------

delay_strings = df1.columns.values[1:]
td = delay_strings.astype(float)
OD_exp = []

for delay in delay_strings:
    OD = df1[delay]
    OD_roi = OD[np.all([(w_ev > roi[0]), (w_ev < roi[1])], axis=0)]
    OD_exp.append(OD)

    # restrict fit parameters to resonable values
    p_init[:-1:3] = abs(p_init[:-1:3])
    p_init[2:-1:3] = abs(p_init[2:-1:3])    

    popt, pcov = curve_fit(hf.fit_lineshapes, w_roi, OD_roi, 
                             p_init, maxfev=1000000, bounds=bounds)

    params.append(popt)
    errors.append(np.sqrt(np.abs(np.diag(pcov))))
    
    # reuse fit result in next iteration
    # p_init = popt

params = np.array(params)
errors = np.array(errors)

#-----------------------------------------------------------------------
#                   Plots
#-----------------------------------------------------------------------

# plt.figure(1)
fig, ((ax1,ax2,ax3),(ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3)

ax1.plot(td, params[:,0])
ax1.set_title('Line strength')
ax1.set_ylabel('T1')

ax2.plot(td, params[:,1])
ax2.set_title('Phase (rad)')

ax3.plot(td, params[:,2])
ax3.set_title('Line width (eV)')

ax4.plot(td, params[:,3])
ax4.set_ylabel('T2')

ax5.plot(td, params[:,4])

ax6.plot(td, params[:,5])

ax7.plot(td, params[:,6])
ax7.set_ylabel('T3')

ax8.plot(td, params[:,7])
ax8.set_xlabel('Time delay (fs)')

ax9.plot(td, params[:,8])

plt.figure(2)
fig, ((ax1,ax2,ax3)) = plt.subplots(1, 3)

ax1.plot(td, params[:,0])
ax2.plot(td, params[:,1])
ax3.plot(td, params[:,2])

ax1.set_title('Line strength')
ax2.set_title('Phase (rad)')
ax3.set_title('Line width (eV)')
ax2.set_xlabel('Time delay (fs)')

ax1.set_xlim([-5,6])
ax1.set_ylim([-0.01,0.22])
ax2.set_xlim([-5,6])
ax2.set_ylim([-2.25,1])
ax3.set_xlim([-5,6])
ax3.set_ylim([0.09,0.33])

plt.show()
#-----------------------------------------------------------------------
#                   Image
#-----------------------------------------------------------------------

# paramdict = {'cmap':'jet', 'shading':'nearest'}
# fontdict = {'color':'w', 'fontsize':10}

# OD_exp = np.array(OD_exp)
# OD_fit = np.zeros((params.shape[0], w_ev.size))
# for p, popt in enumerate(params):
#     OD_fit[p] = hf.fit_lineshapes(w_ev, *popt)

# plt.figure(6)
# fig, ax = plt.subplots(nrows=1,ncols=2, num=5)
# fig.subplots_adjust(right=0.9, left=0.1, top=0.9, bottom=0.15, wspace=0.2)

# im = ax[0].pcolor(w_ev, td, OD_exp, **paramdict, vmin=-0.04, vmax=0.25)
# vmin, vmax = im.get_clim()
# im2 = ax[1].pcolor(w_ev, td, OD_fit, **paramdict, vmin=vmin, vmax=vmax)

# cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
# cbar = fig.colorbar(im2, cax=cbar_ax)
# cbar_ax.set_title('OD')
# cbar.set_ticks([0,0.1,0.2])

# ax[0].set_xlim(photonenergy[[0,-1]])
# ax[1].set_xlim(photonenergy[[0,-1]])
# ax[0].set_ylabel(r'Time delay (fs)')
# ax[0].set_xlabel(r'Photon energy (eV)')
# ax[1].set_xlabel(r'Photon energy (eV)')

# ax[0].text(55.35,-5, r'$T_1$', **fontdict)
# ax[0].text(56.0 ,-5, r'$T_2$', **fontdict)
# ax[0].text(57.32,-5, r'$T_3$', **fontdict)

# ax[0].set_title('Simulation', fontsize=12, pad=1)    
# ax[1].set_title('Fit', fontsize=12, pad=1)  

# plt.show()

#-----------------------------------------------------------------------
#                   Save Fitting Parameters
#-----------------------------------------------------------------------

ndf = pd.DataFrame()
ndf['Time Delays'] = td 
ndf['Line Strength (T1)'] = params[:,0]
ndf['Phase (T1)'] = params[:,1]
ndf['Line Width (T1)'] = params[:,2]

ndf['Line Strength (T2)'] = params[:,3]
ndf['Phase (T2)'] = params[:,4]
ndf['Line Width (T2)'] = params[:,5]

ndf['Line Strength (T3)'] = params[:,6]
ndf['Phase (T3)'] = params[:,7]
ndf['Line Width (T3)'] = params[:,8]

ndf.to_csv(f'fit_params{intensity}.csv', index=False)