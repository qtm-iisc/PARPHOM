"""
    Author: Shinjan Mandal
    Contact: mshinjan@iisc.ac.in
    Github: ShinjanM
    
    Generation of Force constants by finite difference method.
    
    Usage
    -----

    mpiexec -np <no. of cores> get_force_constant.py -i <lammps input file> 
        
        [ Optional arguments: -o Output file, 
                              -d Displacement of each atom,
                              -s turn on symmetry,
                              -t Number of symmetry operations
                              -r Replicate]
    
    Theory
    ------

    We displace every atom by a finite distance and compute the derivative of the 
    Forces on the atoms due to the displacement.

        $\phi_{a,b}^{i,j} = - \frac{dF_{b,j}}{dr_{a,i}}$

"""



from mpi4py import MPI
import numpy as np
from lammps import lammps
import h5py
import argparse
from time import time
import sys
import os

def parseOptions(comm):
    parser = argparse.ArgumentParser(description=print_logo(comm))

    parser.add_argument("-i","--input", help='Input File', type=str)
    parser.add_argument("-o","--output", help='Output File',
                        nargs=argparse.OPTIONAL,
                        default='FORCE_CONSTANTS')
    parser.add_argument("-d","--displacement",help="Displacement of each atom",
                        nargs=argparse.OPTIONAL,default=0.0001,type=float)
    parser.add_argument("-s","--symmetrize",help="Symmetrize the Force constant",
                        default=False,action="store_true")
    parser.add_argument("-t","--iterations",help="Number of iterations to symmetrize",
                        nargs=argparse.OPTIONAL,default=20,type=int)
    parser.add_argument("-r","--replicate",help="Replicate the cell in x, y and z directions. Sample input: '-r 2 3 5'  replicates the cell to create a (2,3,5) supercell. Default is 1 1 1'",
                        nargs=3,type=int,default=[1,1,1])
    args = None
    try:
        if comm.Get_rank() == 0:
            args = parser.parse_args()
    finally:
        args = comm.bcast(args, root=0)

    if args is None:
        exit(0)
    return args



def printProgressBar(i,m,postText,n_bar=20):
    """
        Prints progress of the calculations
    """
    j=i/m
    sys.stdout.write('\r')
    sys.stdout.write(f" {'█' * int(n_bar * j):{n_bar}s} {int(100 * j)}% {postText}")
    sys.stdout.flush()
    
    return


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
                    [0,   0,  zhi-zlo]])
    type_mass = lmp.extract_atom("mass", 2)
    type_ = lmp.gather_atoms("type", 0, 1)
    masses = np.array([type_mass[type_[i]] for i in range(na)], dtype=float)
    at_types = np.array([type_[i] for i in range(na)],dtype=int)
    xp = lmp.extract_atom("x", 3)
    positions = np.array([[xp[i][0], xp[i][1], xp[i][2]] for i in range(na)], dtype=float)
    crys = np.array([np.linalg.solve(lat,positions[i]) for i in range(na)])
    return na, lat, type_mass, at_types, masses, positions,crys




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





def get_forces(lmp, id_, reference, natom, d,at_id,alpha):
    ref = np.copy(reference)
    ref[at_id,alpha] += d
    for i in range(natom):
        lmp.command('set atom {} x {} y {} z {}'.format(id_[i]+1,
                                                        ref[i,0],
                                                        ref[i,1],
                                                        ref[i,2]))
    lmp.command('run 0')
    fp = lmp.extract_atom("f",3)
    forces = np.array([[fp[i][0], fp[i][1], fp[i][2]] for i in range(natom)],dtype=np.float64)
    return forces




def get_force_constant(lammps_input_file,displacement,at_id,alpha,replicate,show_log=False):
    
    cmd_list = ['-log', 'none']
    if not show_log:
        cmd_list += ['-echo', 'none', '-screen', 'none']
    lmp = lammps(cmdargs=cmd_list)
    lmp.file(lammps_input_file)
    
    # Replicate the atoms
    # -------------------

    # For standard potentials this command would work

    lmp.command('replicate {} {} {}'.format(replicate[0],replicate[1],replicate[2]))

    # ---------------------------------------------------------------------
    # For potentials that work between different molecular types,
    # you have to define the criteria for defining the molecular types in 
    # a replicated layout
    # Uncomment the following lines if you want are working with the Drip
    # potential in tBLG and want to define layer 1 as molecule 1 and 
    # layer 2 as molecule 2
    #
    # Refer: https://github.com/lammps/lammps/issues/3047
    #
    # Note that if you are working with multi-layered graphene systems and
    # Drip potential, you need to have cutoffs for multiple layers.
    # ---------------------------------------------------------------------

    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #
    # UNCOMMENT THIS SECTION FOR DRIP POTENTIALS
    #
    #layer_cutoff = 16.5  # change this cutoff depending upon your convenience
    #lmp.command('variable z atom z<%f'%layer_cutoff)
    #lmp.command('group layer1 variable z')
    #lmp.command('group layer2 subtract all layer1')
    #lmp.command('set group layer1 mol 1')
    #lmp.command('set group layer2 mol 2')
    # 
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    lmp.command('run 0')

    na = lmp.get_natoms()
   
    fp = lmp.extract_atom("f",3)
    forces = np.array([[fp[i][0], fp[i][1], fp[i][2]] for i in range(na)], dtype=np.float64)
    
    id_ = lmp.extract_atom("id",0)
    id_ = np.array([id_[i]-1 for i in range(na)], dtype=np.int64)
    
    xp = lmp.extract_atom("x", 3)
    reference = np.array([[xp[i][0], xp[i][1], xp[i][2]] for i in range(na)], dtype=np.float64)
    f1 = get_forces(lmp, id_, reference, na, -displacement/2, at_id, alpha)
    f2 = get_forces(lmp, id_, reference, na,  displacement/2, at_id, alpha)
    fc = -(f2-f1)/displacement
    
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
        if comm.rank==0:
            printProgressBar(i+1,no_of_iterations," symmetrization completed")
            #print("Symmetry iteration %6d/%d completed"%(i+1, no_of_iterations),flush=True)
    f.close()
    return

def acoustic_sum_rule(data,local_atoms):
    """
        Impose acoustic sum rule
    """
    for i in local_atoms:
        for alpha in range(3):
            for beta in range(3):
                col_j = data[i,:,alpha,beta]
                col_j -= np.mean(col_j)
                data[i,:,alpha,beta] = col_j
    return


def inversion_symmetry(data,local_atoms,natom):
    """
        Impose inversion symmetry
    """
    for i in local_atoms:
        for j in range(i,natom):
            for alpha in range(3):
                for beta in range(3):
                    mean = (data[i,j,alpha,beta] + data[j,i,beta,alpha])/2
                    data[i,j,alpha,beta] = mean
                    data[j,i,beta,alpha] = mean
    return

def p2s_map(lat,lammps_input_file, replicate,natom,crys):
    """
        Primitive to supercell map in case of force constants generated in supercells
        Returns the indices of the atoms in the replicated lammps file which are 
        equivalent to the atoms in the original unit cell.
    """

    cmd_list = ['-log', 'none', '-echo', 'none', '-screen', 'none']
    lmp = lammps(cmdargs=cmd_list)
    lmp.file(lammps_input_file)
    lmp.command('run 0')
    lmp.command('replicate {} {} {}'.format(replicate[0],replicate[1],replicate[2]))
    lmp.command('run 0')
    na = lmp.get_natoms()
    id_ = lmp.extract_atom("id",0)
    id_ = np.array([id_[i]-1 for i in range(na)], dtype=np.int64)

    xp = lmp.extract_atom("x", 3)
    xp_rep = np.array([[xp[i][0], xp[i][1], xp[i][2]] for i in range(na)], dtype=np.float64)
    xp_crys = np.array([np.linalg.solve(lat,xp_rep[i]) for i in range(na)])
    p2s_map = np.array(xp_crys-np.tile(crys,(np.prod(replicate),1)))
    p2s_map[np.abs(p2s_map)<1e-8] = 0
    na_arr = np.array([i for i in range(na)])
    p_map = np.where(~p2s_map.any(axis=1))[0]
    return p2s_map, p_map




def print_logo(comm):
    r = comm.Get_rank()
    if r == 0:
        print(" ")
        print(" ")
        print("\t██████╗░░█████╗░██████╗░░█████╗░██████╗░██╗░░██╗░█████╗░░██████╗",flush=True)
        print("\t██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔══██╗██║░░██║██╔══██╗██╔════╝",flush=True)
        print("\t██████╔╝███████║██████╔╝███████║██████╔╝███████║██║░░██║╚█████╗░",flush=True)
        print("\t██╔═══╝░██╔══██║██╔══██╗██╔══██║██╔═══╝░██╔══██║██║░░██║░╚═══██╗",flush=True)
        print("\t██║░░░░░██║░░██║██║░░██║██║░░██║██║░░░░░██║░░██║╚█████╔╝██████╔╝",flush=True)
        print("\t╚═╝░░░░░╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░",flush=True)
        print(" ", flush=True)
        print(" ", flush=True)
        print("\t\t\tForce Constant Generator v1.0",flush=True)
        print(" ", flush=True)
        print(" ", flush=True)
        print("\t\tS. Mandal, I. Maity, H. R. Krishnamurthy, M. Jain",flush=True)
        print(" ", flush=True)
        print(" ", flush=True)
    return





def main():

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    st = time()
    args = parseOptions(comm)
    args.replicate = np.array(args.replicate)
    
    comm.Barrier()

    lammps_input_file = args.input
    na, lat, type_mass, at_types, masses, positions, crys = get_structure_from_lammps(lammps_input_file)
    na_replicated = na*np.prod(args.replicate)


    f = h5py.File(args.output,'w',driver='mpio',comm=comm)

    if np.allclose(args.replicate,np.array([1,1,1]))!=True:
        p1,p2 = p2s_map(lat,lammps_input_file, args.replicate, na,crys)
        dset_p2smap = f.create_dataset('p2s_map',data=p1,dtype=int)
        dset_p_indices = f.create_dataset('p_indices',data=p2,dtype=int)

    dset = f.create_dataset('force_constants',((na_replicated,na_replicated,3,3)),dtype=np.float64)
    
    at_id = [i for i in range(na_replicated)]
    
    local_atoms = distribute(at_id,rank,size)
    
    for j in local_atoms:
        for beta in range(3):
            fc = get_force_constant(lammps_input_file,args.displacement,j,beta,args.replicate)
            dset[j,:,beta] = np.array(fc,dtype='f')
        if rank==0:
            printProgressBar(j+1,len(local_atoms),"force constant generation completed")
    f.close()
    if rank==0: 
        print("\n\nForce constant generated and written in %.3f sec\n"%(time()-st), flush=True)
    comm.Barrier()
    st = time()
    if args.symmetrize:
        symmetrize_fc(args.output,comm,local_atoms,na_replicated,args.iterations)
        if rank==0:
            print("\n\nSymmetrization finished in %.3f sec\n"%(time()-st),flush=True)
    return


if __name__ == '__main__':
    main()
