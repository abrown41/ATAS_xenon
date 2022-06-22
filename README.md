This repository contains the data and analysis scripts for the paper 'Core-resonance
line-shape analysis of atoms undergoing strong-field ionization'.

The required data (both experimental and RMT) is in the `data` directory. 

The `bin` directory houses the analysis scripts used to extract ATAS from the RMT 
dipole files (found in `data/rmt`) and perform the fit. The same analysis scripts
can be used on the experimental data found in `data/experimental`. Currently the
analysis scripts only fit the first peak (T1) of the experimental data.

The `archer_files` directory contains a utility for setting up
calculations, a utility for generating the gaussian profile pulses, and the 
necessary input.conf file to make that work. Also in `archer_files` are utilities 
for extracting the groundstate population at the end of the RMT calculation for each
IR intensity (a single time delay is chosen because the population agrees to \~6s.f. 
across the time delays).

*For users of the RMT suite: a csv file in the correct format can be generated from 
rmt output files using the `bin/makedf.py` utility.
