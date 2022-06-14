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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import helper_functions as hf
import warnings

warnings.filterwarnings("ignore")

args = hf.read_command_line()
energy_shift = args['energy_shift']
intensity = args["IR_intensity"]
literature_linewidth = 0.122  # eV linewidth from literature (Anderson 2001)


def AugerDecayFactor(time_in_au, t_zero=904.1058):
    """
    Prepare the exponential decay factor to simulate Auger decay in the
    simulated dipole. The decay factor is exp(-t/lifetime) where liftime is the
    life time of the resonance (computed from the linewidth

    Parameters
    ----------
    time_in_au : list-like
        array of times at which the dipole has been evaluated in atomic units.
    t_zero : float
        the time at which to start the exponential decay which is the end of
        the XUV pulse: 1 FWHM after the peak of the XUV.
    """
    time = time_in_au - t_zero
    lifetime = 1/(hf.ev_to_au(literature_linewidth)/2)
    decay = np.exp(- time/lifetime)
    decay[time < 0] = 1
    return decay


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
    OD_sim : pd.DataFrame
        dataframe containing the optical density computed at each time delay.
    """
    from scipy.constants import physical_constants
    alpha = physical_constants['fine-structure constant'][0]
    roi = args['roi']

    df1 = pd.read_csv(f"dipole{intensity}.csv")
    df2 = pd.read_csv(f"field{intensity}.csv")


# prepad signal with zeros to improve resolution
    file_length = len(df1)
    desired_length = 2**21
    num_zeros_to_pad = desired_length - file_length

    decay = AugerDecayFactor(df1['Time'])

    dt = df1['Time'][1] - df1['Time'][0]
    w = np.fft.fftfreq(desired_length, dt/(2*np.pi))
    w_ev = hf.au_to_ev(w)
    w_roi = w_ev[np.all([(w_ev > roi[0]), (w_ev < roi[1])], axis=0)]

    OD_sim = {}  # for performance: build dict and then convert to dataframe

    for delay in df1.columns.values[1:]:
        efield = df2[delay]
        xuv_dipole = decay * df1[delay]

        # prepad to improve FT resolution
        field = np.pad(efield, (num_zeros_to_pad, 0),
                       'constant', constant_values=(0, 0))
        dipole = np.pad(xuv_dipole, (num_zeros_to_pad, 0),
                        'constant', constant_values=(0, 0))

        dw = np.fft.fft(dipole)
        ew = np.fft.fft(field)
        OD = -4*np.pi*alpha*w*np.imag(dw/ew)

        OD_roi = OD[np.all([(w_ev > roi[0]), (w_ev < roi[1])], axis=0)]
        OD_sim[delay] = OD_roi

    OD_sim = pd.DataFrame(OD_sim)
    OD_sim['Energy'] = w_roi + energy_shift
    OD_sim.set_index("Energy", inplace=True)

    return OD_sim


def readcsv(fname):
    """
    read the data from csv file 'fname', index it correctly based on the
    filename format
    """
    data = pd.read_csv(fname)
    if fname.startswith('OD'):
        data.set_index('Energy', inplace=True)
    return data


def fitOD(OD_sim):
    """
    Given an energy axis and the optical density for a range of time delays,
    fit the absorption line with the dipole control model, and return the fit
    parameters.

    Parameters
    ----------
    OD_sim : pd.DataFrame
        dataframe containing the optical density computed at each time delay.

    Returns
    -------
    params : pd.DataFrame
        Fit parameters and fitting errors for each time delay
    """
    from scipy.optimize import curve_fit

# fit parameters:[line_strength, phase, line_width, constant_background_term)
#              [  z0,   ϕ,     Γ,             c]
    lower_bounds = [1e-6,  -np.pi, 0.5*literature_linewidth, -1]
    upper_bounds = [np.inf, np.pi, 2*literature_linewidth, 1]
    bounds = (lower_bounds, upper_bounds)

# Initial guess for fitting parameters
#                [z0, ϕ,   Γ,       c]
    p_init = np.array([1, 0, literature_linewidth, 0])
    params = []
    errors = []

    for delay in OD_sim.columns:
        # restrict fit parameters to resonable values
        p_init[:-1:3] = abs(p_init[:-1:3])
        p_init[2:-1:3] = abs(p_init[2:-1:3])

        popt, pcov = curve_fit(hf.fit_lineshapes, OD_sim.index.values,
                               OD_sim[delay], p_init, maxfev=1000000,
                               bounds=bounds)

        params.append(popt)
        errors.append(np.sqrt(np.abs(np.diag(pcov))))

    params = np.array(params)
    errors = np.array(errors)
    params_df = pd.DataFrame()
    params_df['Time Delays'] = [float(time) for time in OD_sim.columns]
    params_df['Line Strength'] = params[:, 0]
    params_df['Phase'] = params[:, 1]
    params_df['Line Width'] = params[:, 2]
    params_df['Background'] = params[:, 3]
    params_df['Line Strength Error'] = errors[:, 0]
    params_df['Phase Error'] = errors[:, 1]
    params_df['Line Width Error'] = errors[:, 2]

    params_df.sort_values(by=['Time Delays'], inplace=True)
    return params_df


def plotParams(params):
    """
    given the time-delays and the fit parameters, plot the line strength, phase
    and line width as a function of time delay

    Parameters
    ----------
    params : pd.DataFrame
        Fit parameters and fitting errors for each time delay

    Returns
    -------
    fig : figure handle for plot
    """

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    axs = {
        'Line Strength': ax1,
        'Phase': ax2,
        'Line Width': ax3
    }
    units = {
        'Line Strength': '',
        'Phase': '(rad)',
        'Line Width': '(eV)'
    }
    for parameter in axs:
        axs[parameter].plot(params['Time Delays'], params[parameter])
        axs[parameter].fill_between(
            params['Time Delays'],
            params[parameter] + params[f'{parameter} Error'],
            params[parameter] - params[f'{parameter} Error'],
            facecolor="gray", alpha=0.35)
        axs[parameter].set_title(f'{parameter} {units[parameter]}')

    ax2.set_xlabel('Time delay (fs)')

    return fig


def getODfit(OD_sim, params):
    """
    Use the fit parameters to reconstruct the OD as a function of energy

    Parameters
    ----------
    OD_sim : pd.DataFrame
        dataframe containing the optical density computed at each time delay.

    params : pd.DataFrame
        Fit parameters and fitting errors for each time delay

    Returns
    -------
    OD_fit : pd.DataFrame
        dataframe containing the reconstructed optical density as a function of
        energy for each time-delay
    """
    OD_fit = {}
    for col in OD_sim.columns:
        row = params.loc[params['Time Delays'] == float(col)]
        popt = [row['Line Strength'].values, row['Phase'],
                row['Line Width'], row['Background']]
        popt = np.reshape(popt, (4,))
        OD_fit[col] = hf.fit_lineshapes(OD_sim.index.values, *popt)

    OD_fit = pd.DataFrame(OD_fit, index=OD_sim.index)
    return(OD_fit)


def plotOD(OD_sim, OD_fit):
    """
    provided with the energy axis, time delay axis and the optical density as a
    function of energy for each time delay, produce a surface plot for both the
    simulated and reconstructed OD.

    Parameters
    ----------
    OD_sim : pd.DataFrame
        dataframe containing the optical density computed at each time delay.
    OD_fit : pd.DataFrame
        dataframe containing the reconstructed optical density as a function of

    Returns
    -------
    fig : figure handle
        figure handle for OD plot
    """
    paramdict = {'cmap': 'jet', 'shading': 'nearest'}

    fig, ax = plt.subplots(nrows=1, ncols=2, num=5)
    fig.subplots_adjust(right=0.9, left=0.1, top=0.9, bottom=0.15, wspace=0.2)

    im = ax[0].pcolor(OD_sim.index.values, [float(x) for x in OD_sim.columns],
                      OD_sim.transpose(), **paramdict, vmin=-0.04, vmax=0.25)
    vmin, vmax = im.get_clim()
    im2 = ax[1].pcolor(OD_sim.index.values, [float(x) for x in OD_fit.columns],
                       OD_fit.transpose(), **paramdict, vmin=vmin, vmax=vmax)

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


def outputData(OD_fit, OD_sim, params):
    """
    Write the simulated and fitted optical densities to file, as well as the
    fitting parameters.

    Parameters
    ----------
    OD_sim : pd.DataFrame
        dataframe containing the optical density computed at each time delay.
    OD_fit : pd.DataFrame
        dataframe containing the reconstructed optical density as a function of
    params : pd.DataFrame
        Fit parameters and fitting errors for each time delay
    """
    OD_sim.to_csv(f'OD{intensity}.csv', index=True)

    OD_fit.to_csv(f'OD_fit{intensity}.csv', index=True)

    params.to_csv(f'fit_params{intensity}.csv', index=False)


def calc_or_read(fname, calcRoutine, **kwargs):
    """ Decide whether to recalculate a particular data structure using
    calcRoutine, or read the pre-calculated data from file, based on user
    input.
    """
    from os.path import isfile

    read_from_file = "n"

    if isfile(fname):
        if args['read_all']:
            read_from_file = 'y'
        else:
            read_from_file = input(f"File {fname} already exists: do you want \
to read data from file? y/n: ")

    if read_from_file == "n":
        data = calcRoutine(**kwargs)
    else:
        data = readcsv(fname)
    return data


def getAllData(intensity):
    """compute, or read from file, the optical density, the fit parameters and
    the fitted optical density"""

    OD_sim = calc_or_read(f"OD{intensity}.csv", calcRoutine=getOD,
                          intensity=intensity)

    params = calc_or_read(f"fit_params{intensity}.csv", calcRoutine=fitOD,
                          OD_sim=OD_sim)

    OD_fit = calc_or_read(f"OD_fit{intensity}.csv", calcRoutine=getODfit,
                          OD_sim=OD_sim, params=params)

    return OD_sim, OD_fit, params


OD_sim, OD_fit, params = getAllData(intensity)

if args['plot']:
    plotParams(params)
    plotOD(OD_sim, OD_fit)

    plt.show()

if args['output']:
    outputData(OD_fit, OD_sim, params)
