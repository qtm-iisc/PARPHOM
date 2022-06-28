import h5py
import matplotlib.pyplot as plt
import numpy as np



#qfile = 'qpt_bands.dat'
qfile = 'q_points_band.dat'

with open(qfile,'r') as f:
    for nodes in f:
        pass

f.close()
nodes = nodes.split()
nodes = np.array([int(i) for i in nodes])
label = [r'$\mathbf{\Gamma}$', r'$\mathbf{M}$', r'$\mathbf{K}$', r'$\mathbf{\Gamma}$']


#f = h5py.File('phbands_21.8_fine_grid.hdf5','r')
f = h5py.File('phbands_21.8_tbwse2_fine_grid.hdf5','r')
g = list(f.keys())

eigval, vel = [], []

for group in g:
    vel.append(f[group]['vel'][:])
    eigval.append(f[group]['eigenvalues'][:])

vel = np.array(vel)
eigval = np.array(eigval)

print(np.shape(vel))

x = [i for i in range(eigval.shape[0])]
#v = [[np.linalg.norm(vel[i,j]) for j in range(vel.shape[1])] for i in range(vel.shape[0])]
v = np.linalg.norm(vel,axis=2)
#v = np.array(v)
print(np.shape(v))
for i in range(eigval.shape[1]):
    plt.plot(x,eigval[:,i],c='b',alpha=0.4,zorder=1)
    plt.scatter(x,eigval[:,i],c=v[:,i],cmap='Reds',zorder=2)

for n in range(len(nodes)):
    plt.axvline(x=nodes[n],linewidth=1,color='k')
plt.axhline(y=0,linewidth=1,color='k')
plt.ylim(-1,100)
plt.xlim(nodes[0],nodes[len(nodes)-1])
plt.ylabel(r"Energy (cm$^{-1}$)",size=20)
plt.xticks(nodes,label,fontsize=22)
plt.tick_params(axis='both',labelsize=20)
cbar = plt.colorbar()
cbar.set_label('Velocity (m/s)',fontsize=18,rotation=270,labelpad=25)
for t in cbar.ax.get_yticklabels():
     t.set_fontsize(20)
plt.tight_layout()
plt.show()

