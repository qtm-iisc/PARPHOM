"""
    Author: Shinjan Mandal
    Contact: mshinjan@iisc.ac.in
    Github: ShinjanM

    Generation of Force constants by finite difference method.
    
    Usage
    -----

    mpiexec -np <no. of cores> get_force_constant -i <lammps input file> 
        
        [ Optional arguments: -o Output file, -d Displacement of each atom]
    
    Theory
    ------

    We displace every atom by a finite distance and compute the derivative of the 
    Forces on the atoms due to the displacement.

        $\phi_{a,b}^{i,j} = -\frac{dF_{b,j}}{dr_{a,i}}$

"""



from mpi4py import MPI
import numpy as np
from lammps import lammps
import h5py
import argparse


def parseOptions(comm):
    parser = argparse.ArgumentParser(description='Print some messages.')

    parser.add_argument("-i","--input", help='Input File', type=str)
    parser.add_argument("-o","--output", help='Output File',
                        nargs=argparse.OPTIONAL,
                        default='FORCE_CONSTANTS')
    parser.add_argument("-d","--displacement",help="Displacement of each atom",
                        nargs=argparse.OPTIONAL,default=0.0001,type=float)
    parser.add_argument("-c","--compression",help="Compress output file",
                        nargs=argparse.OPTIONAL,default=False,type=bool)
    parser.add_argument("-s","--symmetrize",help="Symmetrize the Force constant",
                        nargs=argparse.OPTIONAL,default=False,type=bool)
    args = None
    try:
        if comm.Get_rank() == 0:
            args = parser.parse_args()
    finally:
        args = comm.bcast(args, root=0)

    if args is None:
        exit(0)
    return args


def get_structure_from_lammps(lammps_input_file,show_log=False):

    cmd_list = ['-log', 'none']
    if not show_log:
        cmd_list += ['-echo', 'none', '-screen', 'none']
    lmp = lammps(cmdargs=cmd_list)
    lmp.file(lammps_input_file)
    lmp.command('run 0')
    na = lmp.get_natoms()
    try:
        xlo =lmp.extract_global("boxxlo")
        xhi =lmp.extract_global("boxxhi")
        ylo =lmp.extract_global("boxylo")
        yhi =lmp.extract_global("boxyhi")
        zlo =lmp.extract_global("boxzlo")
        zhi =lmp.extract_global("boxzhi")
        xy =lmp.extract_global("xy")
        yz =lmp.extract_global("yz")
        xz =lmp.extract_global("xz")
    except TypeError:
        xlo =lmp.extract_global("boxxlo", 1)
        xhi =lmp.extract_global("boxxhi", 1)
        ylo =lmp.extract_global("boxylo", 1)
        yhi =lmp.extract_global("boxyhi", 1)
        zlo =lmp.extract_global("boxzlo", 1)
        zhi =lmp.extract_global("boxzhi", 1)
        xy =lmp.extract_global("xy", 1)
        yz =lmp.extract_global("yz", 1)
        xz =lmp.extract_global("xz", 1)
    except UnboundLocalError:
        boxlo, boxhi, xy, yz, xz, periodicity, box_change = lmp.extract_box()
        xlo, ylo, zlo = boxlo
        xhi, yhi, zhi = boxhi
    lat = np.array([[xhi-xlo, xy,  xz],
                    [0,  yhi-ylo,  yz],
                    [0,   0,  zhi-zlo]]).T
    type_mass = lmp.extract_atom("mass", 2)
    type_ = lmp.gather_atoms("type", 0, 1)
    masses = np.array([type_mass[type_[i]] for i in range(na)], dtype=float)
    at_types = np.array([type_[i] for i in range(na)],dtype=int)
    xp = lmp.extract_atom("x", 3)
    positions = np.array([[xp[i][0], xp[i][1], xp[i][2]] for i in range(na)], dtype=float)
    return na, lat, type_mass, at_types, masses, positions


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
    from itertools import islice
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


def get_force_constant(lammps_input_file,displacement,at_id,alpha):
    
    cmd_list = ['-log', 'none', '-echo', 'none', '-screen', 'none']
    lmp = lammps(cmdargs=cmd_list)
    lmp.file(lammps_input_file)
    lmp.command('run 0')
    na = lmp.get_natoms()
   
    fp = lmp.extract_atom("f",3)
    forces = np.array([[fp[i][0], fp[i][1], fp[i][2]] for i in range(na)], dtype=np.float64)
    
    xp = lmp.extract_atom("x", 3)
    reference = np.array([[xp[i][0], xp[i][1], xp[i][2]] for i in range(na)], dtype=np.float64)

    reference[at_id,alpha] -= displacement/2
    id_ = lmp.extract_atom("id", 0)
    id_ = np.array([id_[i]-1 for i in range(na)], dtype=np.int64)
    for i in range(na):
        lmp.command('set atom {} x {} y {} z {}'.format(id_[i]+1,
                                                        reference[i,0],
                                                        reference[i,1],
                                                        reference[i,2]))
    lmp.command('run 0')
    fp1 = lmp.extract_atom("f",3)
    forces1 = np.array([[fp1[i][0], fp1[i][1], fp1[i][2]] for i in range(na)], dtype=np.float64)
    # ---------
    reference[at_id,alpha] += displacement
    id_ = lmp.extract_atom("id", 0)
    id_ = np.array([id_[i]-1 for i in range(na)], dtype=np.int64)
    for i in range(na):
        lmp.command('set atom {} x {} y {} z {}'.format(id_[i]+1,
                                                        reference[i,0],
                                                        reference[i,1],
                                                        reference[i,2]))
    lmp.command('run 0')    
    fp2 = lmp.extract_atom("f",3)
    forces2 = np.array([[fp2[i][0], fp2[i][1], fp2[i][2]] for i in range(na)], dtype=np.float64)
    # ---------
    fc = -(forces2-forces1)/displacement
    return fc


def symmetrize_fc(fc_file,comm,local_atoms,natom,no_of_iterations):
    """
        Impose symmetries
    """
    f = h5py.File(fc_file,'r+',driver='mpio',comm=comm)
    data = f['force_constants']

    for i in range(no_of_iterations):
        acoustic_sum_rule(data,local_atoms)
        inversion_symmetry(data,local_atoms,natom)
    f.close()
    return

def acoustic_sum_rule(data,local_atoms):
    """
        Impose acoustic sum rule
    """
    for i in local_atoms:
        for alpha in range(3):
            for beta in range(3):
                data[i,:,alpha,beta] -= np.mean(data[i,:,alpha,beta])
    return


def inversion_symmetry(data,local_atoms,natom):
    """
        Impose inversion symmetry
    """
    for i in local_atoms:
        for j in range(i+1,natom):
            for alpha in range(3):
                for beta in range(3):
                    mean = (data[i,j,alpha,beta] + data[j,i,beta,alpha])/2
                    data[i,j,alpha,beta] = mean
                    data[j,i,beta,alpha] = mean
    return

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    args = parseOptions(comm)
    lammps_input_file = args.input
    na, lat, type_mass, at_types, asses, positions = get_structure_from_lammps(lammps_input_file)
    at_id = [i for i in range(na)]
    local_atoms = distribute(at_id,rank,size)
    f = h5py.File(args.output,'w',driver='mpio',comm=comm)
    if args.compression:
        dset = f.create_dataset('force_constants',((na,na,3,3)),dtype='f',compression='gzip',compression_opts=9)
    else:
        dset = f.create_dataset('force_constants',((na,na,3,3)),dtype='f')

    cartesian_symb = ['x','y','z']
    for j in local_atoms:
        for beta in range(3):
            fc = get_force_constant(lammps_input_file,args.displacement,j,beta)
            if args.compression:
                with dset.collective:
                    dset[j,:,beta] = np.array(fc,dtype='f')
            else:
                dset[j,:,beta] = np.array(fc,dtype='f')
            print("Rank %d finished displacing atom %d/%d along %s"%(rank, j, local_atoms[len(local_atoms)-1], cartesian_symb[beta]),flush=True)
    f.close()
    if rank==0: 
        print("Force constant generated and written", flush=True)
    comm.Barrier()
    if args.symmetrize:
        symmetrize_fc(args.output,comm,local_atoms,na,20)
    return

if __name__ == '__main__':
    main()
