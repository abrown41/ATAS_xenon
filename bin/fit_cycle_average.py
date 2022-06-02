"""
Reads in fit_params.csv (output from fit_procedure.py) and gs_pop.csv (obtained
from data/pop_GS.* in RMT calculation). Subtracts the cycle-averaged trend from 
the dipole phase to obtain the delay-dependent oscillations. From this, the maximum
oscillation amplitude is obtained for each intensity (neglecting delays <1.75 as in
experiment).

Executed from the directory containing the fit_params.csv files by:
    python fit_cycle_average.py [-p]

The optional argument -p plots three figures:
    1. The delay dependent dipole phase after subtracting the cycye-averaged trend
    2. Maximum amplitude of the oscillating phase at each intensity
    3. Maximum amplitude of the oscillating phase with corresponding Ground State 
       population at end of calculation.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import helper_functions as hf

args=hf.read_command_line()

pop_file=pd.read_csv('gs_pop.csv')
gs_pop = pop_file['Pop_GS']
intensities = pop_file['Intensity']

phase_df = pd.DataFrame()
width_df = pd.DataFrame()
strength_df = pd.DataFrame()
oscil_df =pd.DataFrame()

amplitudes=[]
offset=0

if args["plot"]:
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3)
    ax1.set_xlabel('Time Delay (fs)')
    ax1.set_ylabel('Phase (rad)')
    ax2.set_xlabel('Intensity, $10^{14}$Wcm$^{-2}$')
    ax2.set_ylabel('Phase Oscillation Amplitude (rad)')
    ax3.set_xlabel('Ground State Population')
    ax3.set_ylabel('Phase Oscillation Amplitude (rad)')

for pop,intensity in zip(gs_pop,intensities):
    params = pd.read_csv('fit_params'+str(intensity)+'.csv')
    # print(params['Time Delays'])
    # cut = params[params['Time Delays']==-1.974746].index.values
    cut=[17]
    td = params['Time Delays']
    phase = params['Phase']
    for i in range(len(phase)):
        if phase[i]>1:
            phase[i]-=2*np.pi

    width_df[str(intensity)]=params['Line Width']
    strength_df[str(intensity)]=params['Line Strength']
    phase_df[str(intensity)]=params['Phase']

    avg_phase=hf.smooth_data(phase)
    osscilation = phase-avg_phase
    osscilation = osscilation[cut[0]:]

    oscil_df[str(intensity)]=osscilation

    min_amp = np.amin(osscilation)
    max_amp = np.amax(osscilation)
    del_amp = max_amp-min_amp

    amplitudes.append(del_amp)

    if args["plot"]:
        ax1.plot(td[cut[0]:], osscilation+offset, label=str(intensity))
        ax2.plot(intensity, del_amp, 'o')
        ax3.plot(pop, del_amp, 'o')

oscil_df.insert(0, "Time Delay", td[cut[0]:])
oscil_df.to_csv('extract_oscillations.csv', index=False)

width_df.insert(0, "Time Delay", td)
width_df.to_csv('line_width.csv', index=False)

strength_df.insert(0, "Time Delay", td)
strength_df.to_csv('line_strength.csv', index=False)

phase_df.insert(0, "Time Delay", td)
phase_df.to_csv('phase.csv', index=False)

ndf = pd.DataFrame()
ndf['Intensity'] = intensities
ndf['GS Pop'] = gs_pop
ndf['Amplitude'] = amplitudes
ndf.to_csv('amplitude_gs_pop.csv', index=False)

if args["plot"]:
    plt.show()
