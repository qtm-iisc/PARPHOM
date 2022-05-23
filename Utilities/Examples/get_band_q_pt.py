from pyphutil.pyphutil import moire_phonon_utils
import numpy as np

gamma = np.array([0,0,0])
M = np.array([0.5,0.0,0])
K = np.array([2/3.,1/3.,0])

number_of_points_path_1 = 100
number_of_points_path_2 = 80
number_of_points_path_3 = 120

path = [gamma, number_of_points_path_1, M,
        number_of_points_path_2, K,
        number_of_points_path_3, gamma]

points = moire_phonon_utils()
points.get_qpt_bands(path,output_file='q_points_band.dat')
