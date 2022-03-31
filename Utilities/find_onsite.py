import h5py
import numpy as np

fname = 'bands_0.5_E0-Kpt.hdf5'
nbands = 4 
no_k_pts = 1
index_ = [i for i in range(no_k_pts)]

data = h5py.File(fname, 'r')

eigvals = np.empty((no_k_pts, nbands))

counter = 0
for group in data.keys():
    ds_data = data[group]['eigenvalues']
    eigvals[counter] = ds_data[:]
    counter+=1

print("Offset energy: %0.8f meV"%(1000*(eigvals[0,1]+eigvals[0,2])/2))
