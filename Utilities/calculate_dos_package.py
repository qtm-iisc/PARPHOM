from pymelecutil import pymelecutil as p
import numpy as np

lammps_file = 'lammps_tblg_0.50'
lammps_file_location = '/home/mshinjan/MoireElectronicBands/Examples/tblg_0.50'
natom = 51484
at_types = 1
at_style = 'molecular'
dos_data_file = '/home/mshinjan/Electronic_bands_data/tblg_0.50/DOS/dos_0.5_E0.hdf5'
mesh=np.array([12,12,1])
nbands = 512
en_range = np.array([-0.350,0.350])
en_spacing = 0.0001
full_grid_file = 'full_grid_12X12'

#lammps_file = 'lammps_tblg_21.8'
#lammps_file_location = '/home/mshinjan/MoireElectronicBands/Examples/tblg_21.8'
#natom = 28
#at_types = 1
#at_style = 'molecular'
#dos_data_file = '/home/mshinjan/MoireElectronicBands/Examples/tblg_21.8/dos_21.8_E0_300X300.hdf5'
#mesh=np.array([300,300,1])
#nbands = 28
#en_range = np.array([-12.000,6.500])
#en_spacing = 0.0001
#full_grid_file = 'full_grid_300X300'



moire = p.moire_electron_utils(lammps_file,lammps_file_location,natom,at_types,at_style)
moire.read_lammps()
moire.density_of_states(mesh,dos_data_file,nbands,en_range,en_spacing,
                        method='linear triangulation',full_k_grid_file = full_grid_file,
#                        output_file='dos_21.8_E0.dat')
                        output_file='dos_E0.dat')
