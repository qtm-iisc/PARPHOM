import numpy as np

nqpt = 20
start_k = np.array([0,0,0])
end_k = np.array([0.5,0.,0])

fil = 'q_points.dat'
f = open(fil,'w+')
for i in range(nqpt+1):
    q = start_k + (i)*(end_k-start_k)/(nqpt)
    f.write("%+0.6f  %+0.6f  %+0.6f\n"%(q[0],q[1],q[2]))

