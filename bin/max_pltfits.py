import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser as AP

def read_command_line():
    parser = AP()
    parser.add_argument('-p','--param', type=str,help="parameter for plotting: default is z0 (options phi, gam)",default="z0")
    parser.add_argument('-r','--res', type=str,help="resonance for plotting: default is T1 (options phi, gam)",default="T1")

    return vars(parser.parse_args())


df = pd.read_csv("op.csv")
intensities = [1.3,1.6,1.9,2.2,2.5]

args= read_command_line()
ikey=args["param"]
res=args["res"]
for ii,intens in enumerate(intensities):
#    plt.plot(time_delay_axis[::-1],opt_params[ii,::-1,1],label=str(intens))
    key=res+"_"+ikey+str(intens)
    plt.plot(df["time"],df[key],label=str(intens))
plt.title(ikey + " " + res)
plt.xlabel("time delay (fs)")
plt.legend()
plt.show()
