from pymelecutil import *
import numpy as np

fname = 'lammps.dat'
floc = '/home/mshinjan/Phonons/MoirePhonons/Phonon_Calculations/Examples/tbmos2_9.43'
natom = 42
at_types = 6 
at_style = 'atomic'
mesh = np.array([60,60,1])

utils = moire_electron_utils(file_name=fname,
                              file_location=floc,
                              natom = natom,
                              at_types = at_types,
                              at_style = at_style)

utils.read_lammps()

utils.get_kpt_bz(mesh, use_symmetry=True, GAMMA=True, output_file_name='qpt_dos_9.43.dat',
                 write_full_grid=True, full_grid_file_name="full_grid_60X60")
