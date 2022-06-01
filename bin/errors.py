import matplotlib.pyplot as plt
import pandas as pd


pops = pd.read_csv('pops.txt')

def estimate_numerical_error(pops):
    """
    The error in RMT calculations is manifest in the (lack of) norm
    conservation. Basically, if a calculation is unstable, the norm is very
    different from 1.0. We read these values for each intensity from the
    pops.txt file (which contains output from RMT). The error is then estimated
    from (1-norm)/norm. This is scaled up (somewhat arbitrarily) so that
    error for the lowest intensity (which we would expect to be the most stable
    calculation) is on the 1.5% level.

    We then add a new column to the pops dataframe containing the value of
    this percentage error for each intensity.
    """

    errs =  (pops["worst_pop"] - 1)/pops["worst_pop"]
    scale_factor = 0.015/(errs[0])
    errs *= scale_factor
    pops["Error"] = abs(errs)

    plt.figure()
    plt.plot(pops["Intensity"], 100*pops["Error"])

    plt.xlabel("Intensity, $10^{14}$Wcm$^{-2}$")
    plt.ylabel("Estimated Numerical Error, %")
    plt.show()


estimate_numerical_error(pops)

phase = pd.read_csv('./amplitude_gs_pop.csv')
ebars = pops["Error"] * phase["Amplitude"]
plt.figure()
plt.errorbar(pops["Pop_GS"], phase["Amplitude"], fmt='o', yerr=ebars, capsize=4.0)
plt.xlabel("Ground State Population")
plt.ylabel("Phase Oscillation Amplitude")
plt.show()

plt.figure()
plt.errorbar(pops["Intensity"], phase["Amplitude"], fmt='o', yerr=ebars,
             capsize=4.0)
plt.xlabel("Intensity, $10^{14}$Wcm$^{-2}$")
plt.ylabel("Phase Oscillation Amplitude")
plt.show()

ndf=pd.DataFrame()
ndf['Intensity']=pops["Intensity"]
ndf['Amplitude']=phase['Amplitude']
ndf['GS Pop']=pops["Pop_GS"]
ndf['Error']=pops['Error']
ndf.to_csv('phase_amp_gs_pop.csv', index=False)