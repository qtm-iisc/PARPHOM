from pyphutil.pyphutil import moire_phonon_utils
import numpy as np

fname = 'in.phonon'
floc = '/home/mshinjan/Phonons/MoirePhonons/Force_Constant_Generation/Examples/tblg_21.8/'
mesh = np.array([100,100,1])

nbands = 84
en_range = np.array([-2,1800])
en_spacing = 0.25


data_file = '/home/mshinjan/Phonons/MoirePhonons/Phonon_Calculations/Examples/tblg_21.8/phbands_21.8_tblg_dos.hdf5' 

utils = moire_phonon_utils(lammps_input_file=floc+fname)
utils.read_lammps()
utils.density_of_states(mesh, data_file, nbands, en_range, en_spacing, 
                        method='linear_triangulation',
                        q_grid_file = 'dos_21.8.dat',
                        full_q_grid_file = 'full_grid_100X100' ,
                        output_file='dos.data')
