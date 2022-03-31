import numpy as np

f = 'full_grid_12X12'
mesh = np.array([12,12])
a = np.loadtxt(f, max_rows=mesh[0]*mesh[1])
print(a)
with open(f,'r') as fil:
    lines = fil.readlines()

mapping = lines[mesh[0]*mesh[1]].split()
mapping = np.array([eval(i) for i in mapping]) 
print(mapping)

def get_gridmap(fil,mesh):
    grid = np.loadtxt(f, max_rows=mesh[0]*mesh[1])
    with open(fil,'r') as f:
        lines = f.readlines()
    mapping = lines[mesh[0]*mesh[1]+1].split()
    mapping = np.array([eval(i) for i in mapping])
    return grid, mapping
