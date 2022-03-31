from pymelecutil import moire_electron_utils
import numpy as np

gamma = np.array([0,0,0])
M = np.array([0.5,0.0,0])
K = np.array([2/3.,1/3.,0])

number_of_points_path_1 = 24
number_of_points_path_2 = 16
number_of_points_path_3 = 30

path = [gamma, number_of_points_path_1, M,
        number_of_points_path_2, K,
        number_of_points_path_3, gamma]

points = moire_electron_utils()
points.get_kpt_bands(path,output_file='kpt_bands.dat')
