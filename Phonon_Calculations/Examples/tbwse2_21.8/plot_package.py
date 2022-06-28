from pyphutil.pyphutil import moire_phonon_utils
from pyphutil.pyphutil import plot



lammps_input_file = '/home/mshinjan/Phonons/MoirePhonons/Force_Constant_Generation/Examples/tbwse2_21.8/in.phonon' 
m = moire_phonon_utils(lammps_input_file=lammps_input_file)
m.read_lammps()



p = plot()

label = [r'$\mathbf{\Gamma}$', r'$\mathbf{M}$', r'$\mathbf{K}$', r'$\mathbf{\Gamma}$']
qfile = 'qpt_bands.dat'
datafile = 'phbands_21.8_tbwse2.hdf5'
nqpt = 130
nbands = m.natom*3

p.band_structure(data_file=datafile,
                      nbands = nbands,
                      nqpt = nqpt,
                      qfile=qfile,
                      label = label,
                      save = False)
