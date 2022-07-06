import h5py
import matplotlib.pyplot as plt
import numpy as np

fname = 'phbands_38.2_mos2_KKp.hdf5'
nbands = 42*3
no_k_pts = 181
index_ = [i for i in range(no_k_pts)]

data = h5py.File(fname, 'r')

eigvals = np.empty((no_k_pts, nbands))

counter = 0
for group in data.keys():
    ds_data = data[group]['eigenvalues']
    eigvals[counter] = ds_data[:]
    counter+=1
x = [i for i in range(no_k_pts)]
for i in range(nbands):
    plt.plot(x, eigvals[:,i],c='k')
plt.show()