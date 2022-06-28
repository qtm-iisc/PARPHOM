from pyphutil.pyphutil import moire_phonon_utils

mph = moire_phonon_utils(lammps_input_file='/home/mshinjan/Phonons/MoirePhonons/Force_Constant_Generation/Examples/tbwse2_21.8/in.phonon')

mph.read_lammps()

mph.chirality(ph_data_file='phbands_dos_21.8.hdf5',
              output_file_name='chirality_dos_21.8_tbwse2.hdf5',
              compression=False)
