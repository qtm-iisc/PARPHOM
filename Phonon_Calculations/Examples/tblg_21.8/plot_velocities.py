import h5py
import numpy as np
import matplotlib.pyplot as plt

fil = 'phbands_21.8.hdf5'
f = h5py.File(fil,'r')

groups = list(f.keys())

vel = []
en = []
for g in groups[:20]:
    vel.append(f[g]['vel'][:10]) 
    en.append(f[g]['eigenvalues'][:10])

vel = np.array(vel)
en = np.array(en)
np.set_printoptions(suppress=True, precision=4)
print(en)
