from pyphutil.pyphutil import moire_phonon_utils
import numpy as np

fname = 'lammps.in'
floc = '/home/mshinjan/Phonons/MoirePhonons/Force_Constant_Generation/Examples/tbmos2_38.2/'
mesh = np.array([200,200,1])

utils = moire_phonon_utils(lammps_input_file = floc+fname)
utils.read_lammps()
utils.get_qpt_bz(mesh, 
                 use_symmetry=False, 
                 GAMMA=True, 
                 output_file_name='dos_38.2_200X200.dat',
                 write_full_grid=True, 
                 full_grid_file_name="full_grid_200X200")
