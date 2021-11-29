"""Lynda's script for generating gaussian shaped XUV pulse for input to RMT"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import pandas as pd

def gaussian(x, mu, sig):
	return 1/(np.sqrt(2*np.pi)*sig)*np.exp(-(((x - mu)/sig)**2)/2)

def window(x, mu, sig):
	return np.exp(-((x-mu)**2)/sig)


rmf = np.loadtxt('EField.andrew')
plt.plot(rmf[:,0], rmf[:,3], label='Andrew OG')
# ef = np.amax(rmf[:,3])
# # print(ef)

# t=np.arange(0.005,2500.005,0.01)
# cos = np.sin(2.05*t+0.5)

# mt=(t-483.98)
# o =60
# window = np.exp(-(mt*mt)/o)
# plt.plot(t, ef*window*cos)

time = np.arange(0.005,2500.005,0.01)
sigma = (0.22*41.34)/2.355
centre = (41.34*5)+277.28

# multiply by 41.34 = convert fs to au
# Peak of IR = 277 au
# FWHM = sigma*2.355
# Want FWHM corresponding to 200as => sigma = 0.2/2.355 = 8.27/2.355 a.u.
# 2.05 = Freq. of XUV pulse
# 1.3325 = phase shift to get peaks to match

gaussian1 = gaussian(time, centre, sigma)
gaussian1 *= 0.023860560897644403/np.amax(gaussian1)
time_peak_gaussian = time[np.argmax(gaussian1)]

sine = np.sin(2.05*(time))
pulse = sine*gaussian1

time_peak_pulse = time[np.argmax(pulse)]

pulse = gaussian1*np.sin(2.05*(time+abs(time_peak_pulse-time_peak_gaussian)))

plt.plot(time, pulse, '--', label='Mine') #Pulse Profile

# plt.plot(time, gaussian1, label='Gaussian Envelope')
# gaussian2 = gaussian(time, centre, (0.22*41.34)/2.355)
# gaussian2*= 0.023860560897644403/np.amax(gaussian2)
# plt.plot(time, gaussian2)

plt.legend()
plt.show()
df = pd.DataFrame()

zeros = np.zeros(len(pulse))

df[0] = time
df[1] = zeros
df[2] = zeros
df[3] = pulse


df.to_csv("EField.3", index=False, header=False, sep=' ')
