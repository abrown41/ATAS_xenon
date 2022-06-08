# -*- coding: utf-8 -*-
"""
Uses dipoleX.X.csv and fieldX.X.csv (output from makedf.py) to calculate the OD
for a single IR intensity at each time delay in a time delay scan, perform the
fitting procedure and save the time delay scan, the fits at each time delay and
the fitting parameters to .csv files named ODX.X.csv, OD_fitX.X.csv and
fit_paramsX.X.csv respectively. 'X.X' is the IR intensity.

In calculating the OD, the dipole is tapered with an exponential decay to
simulate Auger Decay.

Executed from the directory containing the dipole and electric field files by:
    python fit_procedure.py -i <IR intensity>
Where the argument -i chooses the IR intensity.

Adapted from a script by Max.
"""

# -----------------------------------------------------------------------
#                    Imports
# -----------------------------------------------------------------------
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as constants
import pandas as pd
import matplotlib.pyplot as plt
import helper_functions as hf

args = hf.read_command_line()
intensity = args["IR_intensity"]

# -----------------------------------------------------------------------
#                   Constants
# -----------------------------------------------------------------------
gamma_xe1 = 0.122  # from literature (Anderson 2001 I think?)
alpha = constants['fine-structure constant'][0]

# set region of interest for fit
roi = slice(*hf.get_roi(hf.get_energy_axis(), erange=(60, 62)))
photonenergy = hf.get_energy_axis()[roi]

# -----------------------------------------------------------------------
#                   Read Data
# -----------------------------------------------------------------------


def getOD(intensity):
    """
    Given an intensity, attempt to read the associated dipole and field files,
    apply the exponential decay scaling to the dipole (to approximate Auger
    decay) and calculate the optical density.

    Parameters
    ----------
    intensity : float or str
        intensity value used in the file name for the dipole?.?.csv and
        field?.?.csv files. 

    Returns
    -------
    w_roi : np.array 
        photon energies in the region of interest (roi)
    OD_sim : pd.DataFrame
        dataframe containing the optical density computed at each time delay
    """

    df1 = pd.read_csv(f"dipole{intensity}.csv")
    df2 = pd.read_csv(f"field{intensity}.csv")

# -----------------------------------------------------------------------
#                   Delay loop
# -----------------------------------------------------------------------

# prepad signal with zeros to improve resolution
    file_length = len(df1)
    slen = 2**21
    nzero = slen - file_length
# prepare the exponential decay factor
# t_zero is 1 FWHM after the peak of the XUV.
    t_zero = 21.87
    time = hf.au_to_fs(df1['Time']) - t_zero
    lifetime = hf.au_to_fs(1/(hf.ev_to_au(gamma_xe1)/2))
    decay = np.exp(- time/lifetime)
    decay[time < 0] = 1
    OD_sim = {}  # for performance: build dict and then convert to dataframe

    for delay in df1.columns.values[1:]:
        efield = df2[delay]
        xuv_dipole = decay * df1[delay]

        # prepad to improve FT resolution
        field = np.pad(efield, (nzero, 0), 'constant', constant_values=(0, 0))
        dipole = np.pad(xuv_dipole, (nzero, 0), 'constant',
                        constant_values=(0, 0))

        # FT to frequency domain
        dw = np.fft.fft(dipole)
        ew = np.fft.fft(field)

        dt = time[1] - time[0]
        w = np.fft.fftfreq(field.shape[0], hf.fs_to_au(dt)/(2*np.pi))
        w_ev = hf.au_to_ev(w)
        w_roi = w_ev[np.all([(w_ev > 59), (w_ev < 62)], axis=0)]

        OD = -4*np.pi*alpha*w*np.imag(dw/ew)
        OD_roi = OD[np.all([(w_ev > 59), (w_ev < 62)], axis=0)]
        OD_sim[delay] = OD_roi

    OD_sim = pd.DataFrame(OD_sim)
    return w_roi, OD_sim


def fitOD(w_roi, OD_sim):
    """
    Given an energy axis and the optical density for a range of time delays,
    fit the absorption line with the dipole control model, and return the fit
    parameters.

    Parameters
    ----------
    w_roi : np.array
        energy axis in the region of interest (roi)
    OD_sim : pd.DataFrame
        dataframe containing the optical density as a function of w_roi for
        each time-delay

    Returns
    -------
    params : np.array
        dimensions (number of time delays, 4)
        the four fit parameters for each time delay.
        params[:, 0] is the line strength
        params[:, 1] is the phase
        params[:, 2] is the line width
        params[:, 3] is the background term
    errors : np.array
        fitting errors corresponding to the params array
    """
# -----------------------------------------------------------------------
#                   Fit setup
# -----------------------------------------------------------------------
# fit parameters:[line_strength, phase, line_width, constant_background_term)
#              [  z0,   ϕ,     Γ,             c]
    lower_bounds = [1e-6,  -np.pi, 0.5*gamma_xe1, -1]
    upper_bounds = [np.inf, np.pi, 2*gamma_xe1, 1]
    bounds = (lower_bounds, upper_bounds)

# Initial guess for fitting parameters
#                [z0, ϕ,   Γ,       c]
    p_init = np.array([1, 0, gamma_xe1, 0])
    params = []
    errors = []

    for delay in OD_sim.columns:
        # restrict fit parameters to resonable values
        p_init[:-1:3] = abs(p_init[:-1:3])
        p_init[2:-1:3] = abs(p_init[2:-1:3])

        popt, pcov = curve_fit(hf.fit_lineshapes, w_roi, OD_sim[delay],
                               p_init, maxfev=1000000, bounds=bounds)

        params.append(popt)
        errors.append(np.sqrt(np.abs(np.diag(pcov))))

        # reuse fit result in next iteration:
        # p_init = popt

    params = np.array(params)
    errors = np.array(errors)
    return params, errors


def plotParams(td, params, errors):
    """
    given the time-delays and the fit parameters, plot the line strength, phase
    and line width as a function of time delay

    Parameters
    ----------
    td : list
        list of time delay values
    params : np.array
        dimensions (number of time delays, 4)
        the four fit parameters for each time delay.
        params[:, 0] is the line strength
        params[:, 1] is the phase
        params[:, 2] is the line width
        params[:, 3] is the background term
    errors : np.array
        fitting errors corresponding to the params array

    Returns
    -------
    fig : figure handle for plot
    """

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(td, params[:, 0])
    ax1.fill_between(td, params[:, 0]+errors[:, 0],
                     params[:, 0]-errors[:, 0],  facecolor="gray", alpha=0.35)
    ax1.set_title('Line strength')

    ax2.plot(td, params[:, 1])
    ax2.fill_between(td, params[:, 1]+errors[:, 1],
                     params[:, 1]-errors[:, 1],  facecolor="gray", alpha=0.35)
    ax2.set_title('Phase (rad)')
    ax2.set_xlabel('Time delay (fs)')

    ax3.plot(td, params[:, 2])
    ax3.fill_between(td, params[:, 2]+errors[:, 2],
                     params[:, 2]-errors[:, 2],  facecolor="gray", alpha=0.35)

    ax3.set_title('Line width (eV)')

    return fig


def getODfit(w_roi, cols, params):
    """
    Use the fit parameters to reconstruct the OD as a function of w_roi

    Parameters
    ----------
    w_roi : np.array
        energy axis in the region of interest (roi)
    cols : list-like
        time delay headings for dataframe
    params : np.array
        dimensions (number of time delays, 4)
        the four fit parameters for each time delay.
        params[:, 0] is the line strength
        params[:, 1] is the phase
        params[:, 2] is the line width
        params[:, 3] is the background term

    Returns
    -------
    OD_fit : pd.DataFrame
        dataframe containing the reconstructed optical density as a function of
        w_roi for each time-delay
    """
    OD_fit = {}
    for col, popt in zip(cols, params):
        OD_fit[col] = hf.fit_lineshapes(w_roi, *popt)

    OD_fit = pd.DataFrame(OD_fit)
    return(OD_fit)


def plotOD(w_roi, td, OD_sim, OD_fit):
    """
    provided with the energy axis, time delay axis and the optical density as a
    function of energy for each time delay, produce a surface plot for both the
    simulated and reconstructed OD.

    Parameters
    ----------
    w_roi : np.array
        energy axis in the region of interest (roi)
    td    : list-like
        time-delay axis
    OD_sim : pd.DataFrame
        dataframe containing the simulated optical density as a function of
        w_roi for each time-delay
    OD_fit : pd.DataFrame
        dataframe containing the reconstructed optical density as a function of
        w_roi for each time-delay

    Returns
    -------
    fig : figure handle
        figure handle for OD plot
    """
    paramdict = {'cmap': 'jet', 'shading': 'nearest'}

    fig, ax = plt.subplots(nrows=1, ncols=2, num=5)
    fig.subplots_adjust(right=0.9, left=0.1, top=0.9, bottom=0.15, wspace=0.2)

    im = ax[0].pcolor(w_roi-5.5, td, OD_sim.transpose(), **paramdict,
                      vmin=-0.04, vmax=0.25)
    vmin, vmax = im.get_clim()
    im2 = ax[1].pcolor(w_roi-5.5, td, OD_fit.transpose(), **paramdict,
                       vmin=vmin, vmax=vmax)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])
    cbar = fig.colorbar(im2, cax=cbar_ax)
    cbar_ax.set_title('OD')
    cbar.set_ticks([0, 0.1, 0.2])

    ax[0].set_ylabel(r'Time delay (fs)')
    ax[0].set_xlabel(r'Photon energy (eV)')
    ax[1].set_xlabel(r'Photon energy (eV)')

    ax[0].set_title('Simulation', fontsize=12, pad=1)
    ax[1].set_title('Fit', fontsize=12, pad=1)

    return(ax)


def outputData(w_roi, td, OD_fit, OD_sim, params, errors, intensity):
    OD_sim.insert(loc=0, column='Energy', value=w_roi)
    OD_sim.to_csv(f'OD{intensity}.csv', index=False)

    OD_fit.insert(loc=0, column='Energy', value=w_roi)
    OD_fit.to_csv(f'OD_fit{intensity}.csv', index=False)

    params_df = pd.DataFrame()
    params_df['Time Delays'] = td
    params_df['Line Strength'] = params[:, 0]
    params_df['Phase'] = params[:, 1]
    params_df['Line Width'] = params[:, 2]
    params_df['Line Strength Error'] = errors[:, 0]
    params_df['Phase Error'] = errors[:, 1]
    params_df['Line Width Error'] = errors[:, 2]

    params_df.to_csv(f'fit_params{intensity}.csv', index=False)


w_roi, OD_sim = getOD(intensity)
params, errors = fitOD(w_roi, OD_sim)
td = [-float(delay) for delay in OD_sim.columns]

OD_fit = getODfit(w_roi, OD_sim.columns, params)

plotParams(td, params, errors)
plotOD(w_roi, td, OD_sim, OD_fit)

plt.show()

outputData(w_roi, td, OD_fit, OD_sim, params, errors, intensity)
