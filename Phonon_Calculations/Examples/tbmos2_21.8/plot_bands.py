import h5py
import matplotlib.pyplot as plt
import numpy as np

<<<<<<< HEAD:Phonon_Calculations/Examples/tbmos2_21.8/plot_bands.py
fname = 'phbands_21.8_mos2.hdf5'
nbands = 42*3
=======
fname = 'phbands_21.8_tblg.hdf5'
nbands = 28*3
>>>>>>> 2c64b51a4f5c46689e30c8cc106bbe247a821212:Phonon_Calculations/Examples/tblg_21.8/plot_bands.py
no_k_pts = 130
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
