from pyphutil.pyphutil import moire_phonon_utils
import numpy as np

fname = 'in.phonon_tersoff'
floc = '../../Force_Constant_Generation/Examples/tblg_21.8/'      # Location of Force Constants
mesh = np.array([15,15,1])

nbands = 84
en_range = np.array([-10,1800])
en_spacing = 0.04169585375


data_file = '../../Phonon_Calculations/Examples/tblg_21.8/phbands_21.8_dos_15X15_tersoff.hdf5' # Phonon Data file with all the q-points

utils = moire_phonon_utils(lammps_input_file=floc+fname)
utils.read_lammps()
utils.density_of_states(mesh, 
                        data_file, 
                        nbands, 
                        en_range, 
                        en_spacing, 
                        method='gaussian',
                        q_grid_file = '../../Phonon_Calculations/Examples/tblg_21.8/dos_21.8_15X15.dat',
                        full_q_grid_file = '../../Phonon_Calculations/Examples/tblg_21.8/full_grid_15X15' ,
                        output_file='dos_15X15_gauss_tersoff.data',
                        func =None,
                        width=5)
