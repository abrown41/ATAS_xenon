"""
Utility for changing the time step and the time delay column labels
of dipoleX.X.csv and fieldX.X.csv (where X.X is the IR intensity).
"""

import pandas as pd

intensities=[1.1, 1.3, 1.5, 1.7, 1.9, 2.2, 2.5]
files=['dipole', 'field']
for intensity in intensities:
	for f in files:
		df = pd.read_csv('int'+str(intensity)+'/'+f+'.csv')
		ndf = pd.DataFrame()
		cols = df.columns.values
		ndf['Time']=df['Time'][::100]
		for col in cols[1:]:
			# delay=float(col)*27.212/41.34
			ndf[col] = df[col][::100]
		ndf.to_csv('../'+f+str(intensity)+'.csv', index=False)