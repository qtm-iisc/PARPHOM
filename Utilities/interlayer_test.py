from pymelecutil import moire_electron_utils as m
import numpy as np

lammps_file_name = 'lammps_tblg_21.8'
lammps_file_location = '/home/mshinjan/MoireElectronicBands/Examples/tblg_21.8'
natom = 28
at_style = 'molecular'
at_types = 1

moire = m(lammps_file_name, lammps_file_location, natom, at_types, at_style)
moire.read_lammps()
box_dim = np.array([[-16,16],
                    [-16,16]])
moire.interlayer_seperation(4, box_dim, vmin=3.48, vmax=3.58)

