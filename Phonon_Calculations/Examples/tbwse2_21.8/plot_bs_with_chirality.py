import h5py
import matplotlib.pyplot as plt
import numpy as np

#plt.rcParams.update({
#  "text.usetex": True,
#  "font.family": "Helvetica"
#})



qfile = 'qpt_bands.dat'

with open(qfile,'r') as f:
    for nodes in f:
        pass

f.close()
nodes = nodes.split()
nodes = np.array([int(i) for i in nodes])
label = [r'$\mathbf{\Gamma}$', r'$\mathbf{M}$', r'$\mathbf{K}$', r'$\mathbf{\Gamma}$']

f = h5py.File('phbands_21.8_tbwse2.hdf5','r')
g = list(f.keys())

eigval, vel = [], []

for group in g:
    vel.append(f[group]['vel'][:])
    eigval.append(f[group]['eigenvalues'][:])

vel.append(f[g[0]]['vel'][:])
eigval.append(f[g[0]]['eigenvalues'][:])


vel = np.array(vel)
eigval = np.array(eigval)

f_chiral = h5py.File('chirality_21.8_tbwse2.hdf5','r')
s_z = f_chiral['chirality'][:,:,2]
s_z = np.append(s_z, [list(s_z[0])], axis=0 )

print(s_z)

print(np.shape(s_z))

x = [i for i in range(eigval.shape[0])]
v = [[np.linalg.norm(vel[i,j]) for j in range(vel.shape[1])] for i in range(vel.shape[0])]
v = np.array(v)
for i in range(eigval.shape[1]):
    plt.plot(x,eigval[:,i],c='gray',alpha=0.4,linewidth=0.6,zorder=1)
    plt.scatter(x,eigval[:,i],c=s_z[:,i],cmap='coolwarm',alpha=(abs(s_z[:,i])/np.max(abs(s_z))),vmax=1,vmin=-1,zorder=2,s=(abs(s_z[:,i])/np.max(abs(s_z)))*50  )
for n in range(len(nodes)):
    plt.axvline(x=nodes[n],linewidth=1,color='k')
plt.axhline(y=0,linewidth=1,color='k')
plt.xlim(nodes[0],nodes[len(nodes)-1])
plt.ylabel(r"Energy (cm$^{-1}$)",size=20)
plt.xticks(nodes,label,fontsize=22)
plt.tick_params(axis='both',labelsize=20)
cbar = plt.colorbar()
cbar.set_label('Chirality',fontsize=18,rotation=270,labelpad=25)
for t in cbar.ax.get_yticklabels():
     t.set_fontsize(20)
plt.show()
