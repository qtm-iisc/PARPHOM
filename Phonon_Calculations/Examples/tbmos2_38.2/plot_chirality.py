import h5py
import matplotlib.pyplot as plt
import numpy as np

#plt.rcParams.update({
#  "text.usetex": True,
#  "font.family": "Helvetica"
#})

import matplotlib



qfile = 'qpt_bands_KKp.dat'

with open(qfile,'r') as f:
    for nodes in f:
        pass

f.close()
nodes = nodes.split()
nodes = np.array([int(i) for i in nodes])
label = [r'$\mathbf{\Gamma}$', r'$\mathbf{M}$', r'$\mathbf{K}$', r'$\mathbf{\Gamma}$', r"$\mathbf{K'}$",r'$\mathbf{M}$']

f = h5py.File('phbands_38.2_mos2_KKp.hdf5','r')
g = list(f.keys())

eigval, vel = [], []

for group in g:
    vel.append(f[group]['vel'][:])
    eigval.append(f[group]['eigenvalues'][:])

vel = np.array(vel)
eigval = np.array(eigval)

f_chiral = h5py.File('chirality_dos_38.2_tbmos2.hdf5','r')
s_z = f_chiral['chirality'][:,:,2]


#figsize_inches = 252/72


labelsize = 40 
fontsize = 35
plt.rc('font', size=35)
matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'


fig,ax = plt.subplots(figsize=(12,12),constrained_layout=True)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)

x = [i for i in range(eigval.shape[0])]
v = [[np.linalg.norm(vel[i,j]) for j in range(vel.shape[1])] for i in range(vel.shape[0])]
v = np.array(v)





for i in range(eigval.shape[1]):
    ax.plot(x,eigval[:,i],c='gray',alpha=0.5,linewidth=0.6,zorder=2)
    im = ax.scatter(x,eigval[:,i],c=s_z[:,i],cmap='RdBu',
                    alpha=(abs(s_z[:,i])/np.max(abs(s_z))),
                    vmax=1,vmin=-1,zorder=3,s=(abs(s_z[:,i])+0.05)*18)



for n in range(len(nodes)):
    ax.axvline(x=nodes[n],linewidth=2.5,color='gray',linestyle='--',zorder=1)

ax.axhline(y=0,linewidth=2.5,color='gray',linestyle='--',zorder=1)
ax.set_xlim(nodes[0],nodes[len(nodes)-1])
ax.set_ylabel(r'$\mathbf{Energy}$ $\mathbf{(cm^{-1})}$', size=labelsize)
ax.set_xticks(nodes,label,fontsize=fontsize)
ax.tick_params(axis='both',which='major',labelsize=labelsize)
cbar = fig.colorbar(im,ax=ax,pad=0.05)
cbar.ax.set_title(r'$\mathbf{s}_{z}^{\mathbf{q}\nu}$',fontsize=fontsize,pad=18)#,rotation=0,labelpad=25)
for t in cbar.ax.get_yticklabels():
     t.set_fontsize(fontsize)
plt.savefig("/home/mshinjan/Final_Images/Phonon_Code/Chiral_MoS2_38.pdf")
#plt.show()

