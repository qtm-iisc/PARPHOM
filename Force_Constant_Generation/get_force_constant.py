# This file is part of PARAPHOM.
# 
# PARAPHOM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PARAPHOM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PARAPHOM.  If not, see <https://www.gnu.org/licenses/>.

"""
@package PARAPHOM.Force_Constant_Generation
@file get_force_constant.py
@brief Generation of force constants by finite difference method using LAMMPS and MPI.

@author Shinjan Mandal, I. Maity, H. R. Krishnamurthy, M. Jain
@copyright GNU General Public License v3.0

This script computes force constants for atomistic systems using the finite difference method. It leverages LAMMPS for force calculations and MPI for parallelization. The code supports cell replication, symmetry operations, and outputs results in HDF5 format.

@par Pre-requisites
- The serial Python wrapper for LAMMPS must be installed. Please follow the official LAMMPS documentation for Python integration: https://docs.lammps.org/Python_module.html

@par Usage
@code{.sh}
mpiexec -np <no. of cores> python get_force_constant.py -i <lammps input file> \
    [-o Output file] [-d Displacement] [-s] [-t Iterations] [-r Replicate]
@endcode

@par Arguments
- -i, --input: Input LAMMPS file (required)
- -o, --output: Output file (default: FORCE_CONSTANTS)
- -d, --displacement: Displacement of each atom (default: 0.0001)
- -s, --symmetrize: Turn on symmetry operations
- -t, --iterations: Number of symmetry iterations (default: 20)
- -r, --replicate: Replicate cell in x, y, z (default: 1 1 1)

@par Theory
The force constant matrix is computed as:
\f[
\phi_{a,b}^{i,j} = - \frac{dF_{b,j}}{dr_{a,i}}
\f]
where atoms are displaced and the resulting forces are used to compute derivatives.

@par Output
- HDF5 file containing the force constant matrix and mapping information.

@par Example
@code{.sh}
mpiexec -np 4 python get_force_constant.py -i in.phonon -o FORCE_CONSTANTS -d 0.0001 -s -t 20 -r 2 2 2
@endcode

@par License
This file is part of PARAPHOM and distributed under the GNU GPL v3.0.
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
    """
    @brief Parse command line options using argparse and broadcast to all MPI ranks.
    @param comm MPI communicator.
    @return argparse.Namespace Parsed arguments.
    """
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
    @brief Prints a progress bar to the console.
    @param i Current iteration.
    @param m Total iterations.
    @param postText Text to display after the progress bar.
    @param n_bar Number of bars in the progress bar.
    """
    j=i/m
    sys.stdout.write('\r')
    sys.stdout.write(f" {'█' * int(n_bar * j):{n_bar}s} {int(100 * j)}% {postText}")
    sys.stdout.flush()
    
    return


def get_structure_from_lammps(lammps_input_file,show_log=False):
    """
    @brief Extract structure information from a LAMMPS input file.
    @param lammps_input_file Path to the LAMMPS input file.
    @param show_log If True, show LAMMPS log output.
    @return Tuple (number of atoms, lattice, type masses, atom types, masses, positions, crystal coordinates)
    """
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
    @brief Distribute a list or array among MPI ranks.
    @param list_to_be_distributed List or numpy array to distribute.
    @param rank Current MPI rank.
    @param size Total number of MPI ranks.
    @return List local to the current rank.
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
    """
    @brief Displace a single atom and compute forces on all atoms.
    @param lmp LAMMPS object.
    @param id_ Atom IDs.
    @param reference Reference positions.
    @param natom Number of atoms.
    @param d Displacement value.
    @param at_id Atom index to displace.
    @param alpha Cartesian direction (0,1,2).
    @return Forces on all atoms after displacement.
    """
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
    """
    @brief Compute force constants by finite difference for a given atom and direction.
    @param lammps_input_file Path to LAMMPS input file.
    @param displacement Displacement value.
    @param at_id Atom index to displace.
    @param alpha Cartesian direction (0,1,2).
    @param replicate Replication factors for supercell.
    @param show_log If True, show LAMMPS log output.
    @return Force constant matrix for the displaced atom and direction.
    """
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
    @brief Impose symmetry operations on the force constant matrix.
    @param fc_file Path to force constant HDF5 file.
    @param comm MPI communicator.
    @param local_atoms List of atom indices local to this rank.
    @param natom Number of atoms.
    @param no_of_iterations Number of symmetry iterations.
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
    @brief Impose the acoustic sum rule on the force constant matrix.
    @param data Force constant matrix.
    @param local_atoms List of atom indices local to this rank.
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
    @brief Impose inversion symmetry on the force constant matrix.
    @param data Force constant matrix.
    @param local_atoms List of atom indices local to this rank.
    @param natom Number of atoms.
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
    @brief Compute the mapping from primitive cell to supercell atom indices.
    @param lat Lattice vectors.
    @param lammps_input_file Path to LAMMPS input file.
    @param replicate Replication factors for supercell.
    @param natom Number of atoms in the primitive cell.
    @param crys Crystal coordinates of atoms in the primitive cell.
    @return Tuple (p2s_map, p_indices)
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
    """
    @brief Print the program logo and credits.
    @param comm MPI communicator.
    """
    r = comm.Get_rank()
    if r == 0:
        print(" ")
        print(" ")
        print("\t██████╗  █████╗ ██████╗ ██████╗ ██╗  ██╗ ██████╗ ███╗   ███╗")
        print("\t██╔══██╗██╔══██╗██╔══██╗██╔══██╗██║  ██║██╔═══██╗████╗ ████║")
        print("\t██████╔╝███████║██████╔╝██████╔╝███████║██║   ██║██╔████╔██║")
        print("\t██╔═══╝ ██╔══██║██╔══██╗██╔═══╝ ██╔══██║██║   ██║██║╚██╔╝██║")
        print("\t██║     ██║  ██║██║  ██║██║     ██║  ██║╚██████╔╝██║ ╚═╝ ██║")
        print("\t╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝")
        print(" ", flush=True)
        print(" ", flush=True)
        print("\t\t\tForce Constant Generator v1.0",flush=True)
        print(" ", flush=True)
        print(" ", flush=True)
        print("\t\tS. Mandal, I. Maity, H. R. Krishnamurthy, M. Jain",flush=True)
        print(" ", flush=True)
        print(" ", flush=True)
        print("\t\t\For any publications arising from this code, please cite:",flush=True)
        print("\t\t\S Mandal, I Maity, H R Krishnamurthy, M Jain, arxiv:2410.21075 ",flush=True)
        print(" ", flush=True)
        print(" ", flush=True)
        print("\t\tThis code is distributed under the GNU General Public License v3.0",flush=True)
        print(" ", flush=True)
    return





def main():
    """
    @brief Main function to generate force constants and optionally symmetrize them.

    Initializes MPI, parses arguments, extracts structure, distributes work, computes force constants,
    writes to HDF5, and applies symmetrization if requested.
    """
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
