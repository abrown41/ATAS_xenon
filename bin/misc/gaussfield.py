"""AB's script for generating gaussian profile XUV pulse to read as input for RMT"""
import numpy as np
import matplotlib.pyplot as plt

rmf = np.loadtxt('rmtfield')
ef = np.max(rmf[:,3])
t=np.arange(0.005,2500.005,0.01)
cos = np.sin(2.05*t+0.5)

mt=(t-483.98)
o =60
window = np.exp(-(mt*mt)/o)

nz=2**18-len(t)
z = np.zeros(nz)

rmf[:,0] = t
rmf[:,3] = ef*window*cos

plt.plot(rmf[:,0], rmf[:,3])
plt.show()
np.savetxt("op.txt",rmf)
