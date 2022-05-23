import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('phbands_21.8_tbwse2_v1.hdf5','r')
g = list(f.keys())

eigval, vel = [], []

for group in g:
    vel.append(f[group]['vel'][:])
    eigval.append(f[group]['eigenvalues'][:])

vel = np.array(vel)
eigval = np.array(eigval)

x = [i for i in range(eigval.shape[0])]
v = [[np.linalg.norm(vel[i,j]) for j in range(vel.shape[1])] for i in range(vel.shape[0])]
v = np.array(v)
for i in range(eigval.shape[1]):
    plt.plot(x,eigval[:,i],c='b',zorder=1)
    plt.scatter(x,eigval[:,i],c=v[:,i],cmap='Reds',zorder=2)
plt.colorbar()
plt.show()
