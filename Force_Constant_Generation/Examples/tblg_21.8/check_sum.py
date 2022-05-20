import h5py 
import numpy as np

f = h5py.File('FORCE_CONSTANTS','r')
fc = f['force_constants'][:]
for i in range(fc.shape[0]):
    for alpha in range(3):
        for beta in range(3):
            print(i,alpha,beta,np.sum(fc[i,:,alpha,beta]))


#for i in range(fc.shape[0]):
#    for j in range(fc.shape[1]):
#        for alpha in range(3):
#            for beta in range(3):
#                print(np.abs(fc[i,j,alpha,beta]-fc[j,i,beta,alpha]))
