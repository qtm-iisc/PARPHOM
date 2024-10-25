from pyphutil.pyphutil import moire_phonon_utils
import numpy as np

fname = 'in.phonon'
floc = '../../Force_Constant_Generation/Examples/tbwse2_2.0/'
mesh = np.array([18,18,1])

nbands = 4902*3
en_range = np.array([-2,1800])
en_spacing = 0.025


data_file = '/home/mshinjan/Twister/examples/Homobilayer_hex/WSe2/Angle_2.00/DOS/phbands_dos_2.00_my_force_with_vel.hdf5'

import h5py
f = h5py.File(data_file,'r')
groups = list(f.keys())
vel = []
for g in groups:
    v = f[g]['vel'][:]
    vel.append((np.square(np.linalg.norm(v,axis=1,ord=2))))

vel = np.array(vel)
print(vel.shape)


utils = moire_phonon_utils(lammps_input_file=floc+fname)
utils.read_lammps()
utils.density_of_states(mesh, 
                        data_file, 
                        nbands, 
                        en_range, 
                        en_spacing, 
                        method='linear_triangulation',
                        q_grid_file = '/home/mshinjan/Twister/examples/Homobilayer_hex/WSe2/Angle_2.00/dos_2.00.dat',
                        full_q_grid_file = 'home/mshinjan/Twister/examples/Homobilayer_hex/WSe2/Angle_2.00/full_grid_18X18' ,
                        output_file='tdos_18X18_twse2.dat',
                        func = vel)
