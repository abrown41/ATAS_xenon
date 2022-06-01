The bin directory houses the analysis scripts used to extract ATAS from the raw
RMT files. The archer_files directory contains a utility for setting up
calculations, a utility for generating the gaussian profile pulses, and the 
necessary input.conf file to make that work. It will need editing to work on your
archer account. Also in the archer_files are utilities for extracting the ground
state population at the end of the calculation for each IR intensity (a single 
time delay is chosen because the population agrees to \~6s.f. across the time 
delays) and the final pop_all which can be used to estimate the error in the
calculation.

The required data (EField and dipole files) are in the data directory. At the moment,
the analysis scripts need tweaked to work with these files.