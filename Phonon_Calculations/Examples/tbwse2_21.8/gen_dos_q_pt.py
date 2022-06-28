from pyphutil.pyphutil import moire_phonon_utils
import numpy as np

fname = 'in.phonon'
floc = '/home/mshinjan/Phonons/MoirePhonons/Force_Constant_Generation/Examples/tbwse2_21.8/'
mesh = np.array([100,100,1])

utils = moire_phonon_utils(lammps_input_file = floc+fname)
utils.read_lammps()
utils.get_qpt_bz(mesh, 
                 use_symmetry=True, 
                 GAMMA=True, 
                 output_file_name='dos_21.8.dat',
                 write_full_grid=True, 
                 full_grid_file_name="full_grid_100X100")
