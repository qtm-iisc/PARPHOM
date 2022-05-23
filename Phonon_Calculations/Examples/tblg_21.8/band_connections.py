import h5py
import numpy as np

f = h5py.File('phbands_21.8_tblg_a.hdf5','r')

groups = list(f.keys())

evec = []

for g in groups:
    evec.append(f[g]['evec'][:])

evec = np.array(evec)

connections = np.zeros((evec.shape[0]+1, evec.shape[1]),dtype=int)
connections[0] = np.array([i for i in range(evec.shape[2])],dtype=int)
connections[len(connections)-1] = connections[0]
print(connections)

for i in range(1,evec.shape[0]):
    e0 = evec[i-1][connections[i-1]]
    e1 = evec[i]
    for j in range(evec.shape[1]):
        connections[i,j] = np.argmax(np.array([abs(np.vdot(e0[k],e1[j])) for k in 
                                        range(evec.shape[2])]))
    print(connections[i])
    print(len(np.unique(connections[i])))
