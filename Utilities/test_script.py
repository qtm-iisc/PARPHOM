import numpy as np
import bz_integration as b

natoms = 52000
nbands = 512

E = np.arange(-250,250,0.1, dtype=np.double)
loc_tri = np.array([[i for i in range(j,j+3)] for j in range(120)])
freq = np.array([np.random.random(natoms) for i in range(nbands)], dtype=np.double)
func = np.array([[1 for i in range(natoms)] for j in range(nbands)], dtype=np.double)

dos, ndos = b.linear_triangulation.sum_over_triangles(nbands, loc_tri, len(E), E, freq, func,0)
