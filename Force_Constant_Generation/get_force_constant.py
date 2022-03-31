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
from iofile import get_structure_from_lammps
from distribute_lists import distribute
import h5py
import argparse


def parseOptions(comm):
    parser = argparse.ArgumentParser(description='Print some messages.')

    parser.add_argument("-i","--input", help='Input File', type=str)
    parser.add_argument("-o","--output", help='Output File',
                        nargs=argparse.OPTIONAL,
                        default='FORCE_CONSTANTS')
    parser.add_argument("-d","--displacement",help="Displacement of each atom",
                        nargs=argparse.OPTIONAL,default=0.01,type=float)
    args = None
    try:
        if comm.Get_rank() == 0:
            args = parser.parse_args()
    finally:
        args = comm.bcast(args, root=0)

    if args is None:
        exit(0)
    return args


def get_force_constant(lammps_input_file,displacement,at_id,alpha):
    
    cmd_list = ['-log', 'none', '-echo', 'none', '-screen', 'none']
    lmp = lammps(cmdargs=cmd_list)
    lmp.file(lammps_input_file)
    lmp.command('run 0')
    na = lmp.get_natoms()
   
    fp = lmp.extract_atom("f",3)
    forces = np.array([[fp[i][0], fp[i][1], fp[i][2]] for i in range(na)], dtype=float)
    
    xp = lmp.extract_atom("x", 3)
    reference = np.array([[xp[i][0], xp[i][1], xp[i][2]] for i in range(na)], dtype=float)

    reference[at_id,alpha] -= displacement/2
    id_ = lmp.extract_atom("id", 0)
    id_ = np.array([id_[i]-1 for i in range(na)], dtype=int)
    for i in range(na):
        lmp.command('set atom {} x {} y {} z {}'.format(id_[i]+1,
                                                        reference[i,0],
                                                        reference[i,1],
                                                        reference[i,2]))
    lmp.command('run 0')
    fp1 = lmp.extract_atom("f",3)
    forces1 = np.array([[fp1[i][0], fp1[i][1], fp1[i][2]] for i in range(na)], dtype=float)
    # ---------
    reference[at_id,alpha] += displacement
    id_ = lmp.extract_atom("id", 0)
    id_ = np.array([id_[i]-1 for i in range(na)], dtype=int)
    for i in range(na):
        lmp.command('set atom {} x {} y {} z {}'.format(id_[i]+1,
                                                        reference[i,0],
                                                        reference[i,1],
                                                        reference[i,2]))
    lmp.command('run 0')    
    fp2 = lmp.extract_atom("f",3)
    forces2 = np.array([[fp2[i][0], fp2[i][1], fp2[i][2]] for i in range(na)], dtype=float)
    # ---------
    fc = -(forces2-forces1)/displacement


    return fc



def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    args = parseOptions(comm)
    lammps_input_file = args.input
    na, lat, type_mass, at_types, asses, positions = get_structure_from_lammps(lammps_input_file)
    at_id = [i for i in range(na)]
    local_atoms = distribute(at_id,rank,size)
    f = h5py.File(args.output,'w',driver='mpio',comm=MPI.COMM_WORLD)
    dset = f.create_dataset('force_constants',((na,na,3,3)),dtype='f')
    for j in local_atoms:
        for alpha in range(3):
            fc = get_force_constant(lammps_input_file,args.displacement,j,alpha)
            dset[j,:,alpha] = fc
            print("Rank %d finished %d:%d"%(rank, j, alpha),flush=True)
    f.close()

if __name__ == '__main__':
    main()
