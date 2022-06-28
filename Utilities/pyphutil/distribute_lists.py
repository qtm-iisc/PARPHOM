"""
    Written by Shinjan Mandal as an utility for Moire Phonons
    Github repo: https://github.com/ShinjanM/MoirePhonons
    Email : mshinjan@iisc.ac.in
"""

from itertools import islice
import numpy as np

def distribute(list_to_be_distributed, rank, size):
    """
        Distributes a list optimally among processors

        Input:
               list_to_be_distributed: list or numpy array to be distributed
               rank: rank of the matrix where the sliced list is to be 
                     distributed to.
               size: total number of processors being used

        Output:
               local_list: list local to the rank
    """
    number_of_local_points = len(list_to_be_distributed)//size
    rem = len(list_to_be_distributed)%size
    local_list = []
    if rank > (size-rem-1):
        number_of_local_points += 1
        local_list = list(islice(list_to_be_distributed, 
                                 rank*number_of_local_points-(size-rem),
                                 (rank+1)*number_of_local_points-(size-rem)))
    else:
        local_list = list(islice(list_to_be_distributed,
                                 rank*number_of_local_points,
                                 (rank+1)*number_of_local_points))
    return local_list



def get_color_and_keys(rank,size,ngroups):
    """
        Get the color and keys for splitting up a MPI communicator into MPI groups

        rank :: rank of the processor under the MPI communicator to be split
        size :: size of the MPI communicator to be split
        ngroups :: number of groups the communicator is to be split into

    """

    neach = size//ngroups
    remeach = size%ngroups
    color_arr = np.array([neach for i in range(ngroups)])
    for j in range(remeach):
        color_arr[ngroups-j-1] += 1
    color_arr = np.cumsum(color_arr)
    color = 0
    key = 0
    while color<ngroups:
        if color_arr[color]>rank:
            key = rank
            break
        else:
            color +=1
    return color, key
