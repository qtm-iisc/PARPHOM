import h5py
import numpy as np

f = h5py.File('phbandsymm.hdf5','r')

groups = list(f.keys())

val = []
evec = []
vel = []
for g in groups:
    val.append(f[g]['eigenvalues'][:])
    evec.append(f[g]['evec'][:])
    vel.append(f[g]['vel'][:])
val = np.array(val)
vel = np.array(vel)
evec = np.array(evec)

S = np.array([[ np.cos(np.pi/3), np.sin(np.pi/3), 0],
              [-np.sin(np.pi/3), np.cos(np.pi/3), 0],
              [               0,               0, 1]])

print(val.shape)

for i in range(val.shape[1]):
    print("%02d   %+02.4f  %+02.4f +%02.4f     %+0.6f"%(i, val[0,i], val[1,i], val[2,i], val[0,i]-val[1,i]))


#for i in range(vel.shape[1]):
#    a = vel[1,i]
#    print("%02d   %+02.3f  %+02.3f  %+02.3f    %+02.3f  %+02.3f  %+02.3f       %0.3f   %0.3f  %0.6f"%(i, vel[0,i,0], vel[0,i,1], vel[0,i,2], a[0], a[1], a[2], np.linalg.norm(vel[0,i]), np.linalg.norm(a), np.linalg.norm(vel[0,i]) - np.linalg.norm(a)) )
