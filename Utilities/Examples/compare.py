import matplotlib.pyplot as plt
import numpy as np


plt.rcParams['text.usetex'] = True

d0k = np.loadtxt('dos.data')#'0Kdos.dat')
plt.fill_between(d0k[:,0],d0k[:,1]*(4.136*2),0,interpolate=True,label=r"$0 K$",alpha=0.6,color='r',linewidth=2,zorder=2)
d300k = np.loadtxt('300Kdos.dat')
plt.fill_between(d300k[:,0],d300k[:,1]*(4*15*15),0,interpolate=True,label=r"$300 K$",color='b',zorder=1)
plt.legend()
plt.xlim(0,1800)
plt.ylim(0,)
plt.xlabel(r"Energy (cm$^{-1}$)",fontsize=20)
plt.ylabel(r"DOS",fontsize=20)
plt.legend(fontsize=20,loc=2,shadow=True,fancybox=True)
plt.tick_params(axis='both',labelsize=20)
plt.tight_layout()
plt.show()
