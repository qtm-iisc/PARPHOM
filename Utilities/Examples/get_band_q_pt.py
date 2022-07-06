from pyphutil.pyphutil import moire_phonon_utils
import numpy as np

gamma = np.array([0,0,0])
M = np.array([0.5,0.0,0])
K = np.array([2/3.,1/3.,0])
Kp = np.array([1/3,2/3,0])
Mp = np.array([0,0.5,0])
number_of_points_path_1 = 40
number_of_points_path_2 = 30
number_of_points_path_3 = 50
number_of_points_path_4 = 50
number_of_points_path_5 = 30

path = [gamma, number_of_points_path_1, M,
        number_of_points_path_2, K,
        number_of_points_path_3, gamma,
        number_of_points_path_4, Kp,
        number_of_points_path_5, Mp]

points = moire_phonon_utils()
points.get_qpt_bands(path,output_file='qpt_bands_KKp.dat')
