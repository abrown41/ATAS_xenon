"""
This script reads both the expec_z_all.* and EField.Xe* files from the individual directories,
and puts them together into two dataframes, which are output to csv files: field.csv and dipole.csv.
These csv files can then be read by the fit_script utility to generate the OD and perform the fitting.
The directories containing the files from the RMT calculation are to be named delay_+_X.XXX or delay_-_X.XXX 
"""
import pandas as pd
import glob
import helper_functions as hf

args = hf.read_command_line()
filelist=args["dirs"]
filelist= hf.order_files(filelist)

ndf_dipole = pd.DataFrame()
ndf_field  = pd.DataFrame()

for f,t in filelist:
    fname=glob.glob(f+"/expec_z_all.*")[0]
    df1 = hf.grab_data(fname)
    ndf_dipole[str(t)] = df1['0001_z']

    fname=glob.glob(f+"/EField.Xe*")[0]
    df2 = hf.grab_data(fname)
    ndf_field[str(t)] = df2['0001_z']

ndf_dipole.insert(0, "Time", df1['Time'])
ndf_dipole.drop(index=1)
ndf_dipole.to_csv('dipole.csv',index=False)

ndf_field.insert(0, "Time", df2['Time'])
ndf_field.drop(index=1)
ndf_field.to_csv('field.csv',index=False)

