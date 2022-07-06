from pyphutil.pyphutil import moire_phonon_utils

mph = moire_phonon_utils(lammps_input_file='/home/mshinjan/Phonons/MoirePhonons/Force_Constant_Generation/Examples/tbmos2_38.2/lammps.in')

mph.read_lammps()

mph.chirality(ph_data_file='phbands_38.2_mos2_KKp.hdf5',
              output_file_name='chirality_dos_38.2_tbmos2.hdf5',
              nbands=42*3)
