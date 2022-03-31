from pymelecutil import *
import numpy as np

fname = 'lammps_tblg_0.95'
floc = '/home/mshinjan/Electronic_bands_data/tblg_0.95'
natom = 14284
at_types = 1 
at_style = 'molecular'
mesh = np.array([40,40,1])

utils = moire_electron_utils(file_name=fname,
                              file_location=floc,
                              natom = natom,
                              at_types = at_types,
                              at_style = at_style)

utils.read_lammps()

utils.get_kpt_bz(mesh, use_symmetry=True, GAMMA=True, output_file_name='dos_0.95.dat',
                 write_full_grid=True, full_grid_file_name="full_grid_40X40")
