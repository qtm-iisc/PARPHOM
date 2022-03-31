"""
    Distributes a list among the available processors optimally
"""

from itertools import islice

def distribute(list_to_be_distributed, rank, size):
    """
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

