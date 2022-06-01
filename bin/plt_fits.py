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
    return vars(parser.parse_args())

args=read_command_line()
labels = ['1.1', '1.3','1.5','1.7', '1.9', '2.2', '2.5']
fig, (ax1,ax2,ax3) = plt.subplots(1, 3)

for i,f in enumerate(args["dirs"]):
    print(f)
    p=Path(f)
    params = pd.read_csv(p / 'fit_params.csv')
    td = params['Time Delays']

    ax1.plot(td, params['Line Strength'], label=labels[i])
    ax1.set_title('Line strength')

    phase=params['Phase']
    for j in range(len(phase)):
        if phase[j]>1:
            phase[j]-=2*np.pi
    ax2.plot(td, phase, label=labels[i])
    ax2.set_title('Phase (rad)')
    ax2.set_xlabel('Time delay (fs)')

    ax3.plot(td, params['Line Width'], label=labels[i])
    ax3.set_title('Line width (eV)')

handles, labels = ax3.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(labels))
plt.show()

