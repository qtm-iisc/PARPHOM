import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('phbands_21.8_tbwse2.hdf5','r')
g = list(f.keys())

eigval, vel = [], []

for group in g:
    vel.append(f[group]['vel'][:])
    eigval.append(f[group]['eigenvalues'][:])

vel = np.array(vel)
eigval = np.array(eigval)

f_chiral = h5py.File('chirality_21.8_tbwse2.hdf5','r')
s_z = f_chiral['chirality'][:,:,2]

print(np.shape(s_z))

print(s_z)

x = [i for i in range(eigval.shape[0])]
v = [[np.linalg.norm(vel[i,j]) for j in range(vel.shape[1])] for i in range(vel.shape[0])]
v = np.array(v)
for i in range(eigval.shape[1]):
    plt.plot(x,eigval[:,i],c='k',alpha=0.5,linewidth=0.8,zorder=1)
    plt.scatter(x,eigval[:,i],c=s_z[:,i],cmap='coolwarm',alpha=(abs(s_z[:,i])/np.max(abs(s_z))),vmax=1,vmin=-1,zorder=2,s=(abs(s_z[:,i])/np.max(abs(s_z)))*10  )
plt.colorbar()
plt.show()
