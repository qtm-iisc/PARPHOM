from pymelecutil import moire_electron_utils as m
import numpy as np

angle = '0.50'

lammps_file_name = 'lammps_tblg_'+angle
lammps_file_location = '/home/mshinjan/Electronic_bands_data/Strain'
natom = 51484
at_style = 'molecular'
at_types = 1

moire = m(lammps_file_name, lammps_file_location, natom, at_types, at_style)
moire.read_lammps()

equib_latcon = 1.42
clim = [-0.1, 0.1]
plot_limits = None #[600,600]
neighbors = 3

heading = r'Angle $%s^{\circ}$'%angle
output = lammps_file_location+'/strain_%s_layer2.png'%angle
moire.strain(equib_latcon, neighbors, heading=heading, plot_limits=plot_limits, save=False, dpi=600, output_file_name=output)

