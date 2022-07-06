"""
    Written by Shinjan Mandal as an utility for Moire Phonons
    Github repo: https://github.com/ShinjanM/MoirePhonons
    Email : mshinjan@iisc.ac.in
"""
import numpy as np
from pyphutil.distribute_lists import distribute
import h5py
import matplotlib.pyplot as plt

class moire_phonon_utils():
    """
        Utility routines for Moire Phonons
    """
    def __init__(self,
                 lammps_input_file=None):
        self.lammps_input_file = lammps_input_file
        self.natom = None
        self.lat = None
        self.at_types = None
        self.mass = None
        self.rec_lat = None
        self.real_pos = None
        self.crys_pos = None
        self.neigh_list = None

    def read_lammps(self, show_log=False):
        """
            Read lammps data file and extract data
        """
        from lammps import lammps
        cmd_list = ['-log', 'none']
        if not show_log:
            cmd_list += ['-echo', 'none', '-screen', 'none']
        lmp = lammps(cmdargs=cmd_list)
        lmp.file(self.lammps_input_file)
        lmp.command('run 0')
        self.natom = lmp.get_natoms()
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
        self.lat = np.array([[xhi-xlo, xy,  xz],
                             [0,  yhi-ylo,  yz],
                             [0,   0,  zhi-zlo]])
        type_mass = lmp.extract_atom("mass", 2)
        self.at_types = lmp.gather_atoms("type",0,1)
        self.mass = np.array([type_mass[self.at_types[i]] for i in range(self.natom)], dtype=float)
        self.at_types = np.array([self.at_types[i] for i in range(self.natom)],dtype=int)
        xp = lmp.extract_atom("x", 3)
        self.real_pos = np.array([[xp[i][0], xp[i][1], xp[i][2]] for i in range(self.natom)], dtype=float)

        self.rec_lat = np.linalg.inv(self.lat.T)*2*np.pi
        
        self.crys_pos = np.array([np.linalg.solve(self.lat,self.real_pos[i])
                                        for i in range(self.natom)])

        return


    def get_qpt_bz(self, 
                   mesh, 
                   use_symmetry=True, 
                   GAMMA=True,
                   output_file_name='kpt_dos.dat',
                   write_full_grid=True,
                   full_grid_file_name='kpt_full.dat'):
        """
            Generates the points on the q-grid mesh you want to sample. 
    
            Input:  Enter the mesh as a numpy array
                    Example: mesh = np.array([5,5,1]), which denotes 
                             a 5X5X1 mesh
                    
                    GAMMA: True for gamma centred mesh, False otherwise  

                    use_symmetry: True/False
            
            Output:
                    output_file: File containing the list of q-points in the requested 
                                 grid in reciprocal lattice.
                                 If symmetry is on, this file stores the reduced grid
                                 information
                    full_grid_file : If symmetry is on, this file stores all the points
                                     in the reciprocal lattice in a n X n grid.
                                     If symmetry is off this file is not created.
        """
        
        import spglib

        cell = (self.lat.T, self.crys_pos, self.at_types)
        if GAMMA==True:
            # Γ centered Mesh
            mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0,0,0])
            if use_symmetry:
                u, inverse = np.unique(mapping, return_inverse=True)
                q_points_crys = (grid[np.unique(u)])/np.array(mesh,dtype=float) 
                full = grid/mesh
                if write_full_grid:
                    f = open(full_grid_file_name,'w')
                    for i in range(len(full)):
                        f.write('%+.6f\t%+.6f\t%+.6f\n'%(full[i,0],
                                                       full[i,1],
                                                       full[i,2]))
                    for i in range(len(mapping)):
                        f.write('%d\t'%mapping[i])
                    f.write('\n')
            else:
                q_points_crys = grid/mesh
        else:
            # Γ shifted Mesh
            mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0.5,0.5,0])
            if use_symmetry:
                u, inverse = np.unique(mapping, return_inverse=True)
                q_points_crys = (grid[np.unique(u)]+[0.5,0.5,0.0])/np.array(mesh,dtype=float)
                full = (grid+np.array([0.5,0.5,0]))/mesh
                if write_full_grid:
                    f = open(full_grid_file_name,'w')
                    for i in range(len(full)):
                        f.write('%+.6f\t%+.6f\t%+.6f\n'%(full[i,0],
                                                       full[i,1],
                                                       full[i,2]))
                    for i in range(len(mapping)):
                        f.write('%d\t'%mapping[i])
                    f.close()
            else:
                q_points_crys = (grid+np.array([0.5,0.5,0.0]))/mesh
        nkpt = mesh[0]*mesh[1]*mesh[2]
        if use_symmetry:
            weights = np.array([0 for i in range(len(q_points_crys))])
            for i in inverse:
                weights[i] += 1
        else:
            weights = np.array([1 for i in range(len(q_points_crys))])
        weights = weights/nkpt
        f = open(output_file_name,'w')
        for i in range(len(q_points_crys)):
            f.write('%+.6f\t%+.6f\t%+.6f\t\t%+.8f\n'%(q_points_crys[i,0],
                                             q_points_crys[i,1],
                                             q_points_crys[i,2],
                                             weights[i]))
        return



    def get_qpt_bands(self, points,output_file='k_points.dat'):
        """
            Generates a list of points along the path specified for plotting the 
            band structure

            Input: 
                  points: list of q-points along which the bands are plotted
            
                  The input style of points is - 
                  [..., q(i), no of points between q(i) 
                                        and q(i+1), q(i+1), ...] 
                  where q(n) is the n-th q-point, in fractional coordinates
                                            
        
            Output: 
                  q_pt: list of the q-point path in fractional coordinates
        """
        q_pt, nodes = [], [0,]
        for i in range(0,len(points)-1,2):
            qi = points[i]
            qiplus1 = points[i+2]
            num = points[i+1]
            nodes.append(num+nodes[len(nodes)-1])
            for i in range(num):
                q_pt.append(qi+((qiplus1-qi)*i/num))
        if (np.allclose(points[0],points[len(points)-1])==False):
            q_pt.append(points[len(points)-1])
        q_pt = np.array(q_pt)
        f = open(output_file,'w')
        for i in range(len(q_pt)):
            f.write('%+.6f\t%+.6f\t%+.6f\n'%(q_pt[i,0],q_pt[i,1],q_pt[i,2]))
        for i in range(len(nodes)):
            f.write('%d\t'%nodes[i])
        return


    def read_full_grid(self, mesh, full_grid_file):
        """
            Read the full k grid.
        """
        grid = np.loadtxt(full_grid_file, max_rows=mesh[0]*mesh[1]*mesh[2])
        with open(full_grid_file,'r') as f:
            lines = f.readlines()
        mapping = lines[mesh[0]*mesh[1]*mesh[2]].split()
        mapping = np.array([eval(i) for i in mapping])
        return grid, mapping


    def density_of_states(self,
                          mesh,
                          data_file,
                          nbands,
                          en_range,
                          en_spacing = 1e-3, #eV
                          method='linear_triangulation',
                          q_grid_file = None,
                          width = 1e-4, #eV
                          full_q_grid_file=None,
                          output_file='dos.dat',
                          func = None):
        """
            Computes the density of states of a system
        """

        
        from bz_integration import bz_integration as bzi
        from mpi4py import MPI
        from pyphutil.progress import printProgressBar

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()  
        size = comm.Get_size()
        E = np.arange(en_range[0],en_range[1],en_spacing)

        if (method=='linear_triangulation'):
            
            from scipy.spatial import Delaunay

            moire_area = abs(np.linalg.det(self.lat[0:2,0:2]))*0.01 # in nm^2
            BZ_area = abs(np.linalg.det(self.rec_lat[0:2,0:2]))*100 # in nm^-2

            grid, mapping = self.read_full_grid(mesh, full_q_grid_file)
            
            Freq = np.empty((len(np.unique(mapping)),nbands))
            
            if func is None:
                func = np.ones(np.shape(Freq),dtype=np.float64)
            elif func.shape!=Freq.shape:
                if rank==0:
                    print("Wrong size for the array func",flush=True)
                exit()

            data = h5py.File(data_file,'r')
            counter = 0
            for group in data.keys():
                ds_data = data[group]['eigenvalues']
                Freq[counter] = ds_data[:]              # in cm^{-1}
                counter += 1
            
            for i in range(len(grid)):
                grid[i] = np.matmul(self.rec_lat.T, grid[i])
            
            grid = np.delete(grid,2,axis=1)
            tri = Delaunay(grid)
            up = list(np.unique(mapping))

            triangles = tri.simplices
            local_triangles_reduced = []
            local_triangles = distribute(triangles,rank,size)
            
            for i in range(len(local_triangles)):
                local_triangles_reduced.append([up.index(mapping[local_triangles[i][0]]), \
                                                up.index(mapping[local_triangles[i][1]]), \
                                                up.index(mapping[local_triangles[i][2]])])
                if rank==0:
                    printProgressBar(i+1,len(local_triangles)," of triangles processed")
                #print("Rank %d appended %d/%d triangles."%(rank,i+1,len(local_triangles)),flush=True)
            comm.Barrier()
            
            if rank==0:
                print("\nStarting DOS calculations",flush=True)

            triangles_reduced = []
            triangles_reduced = comm.gather(local_triangles_reduced, root=0)

            if rank==0:
                triangles_reduced = [i for j in triangles_reduced for i in j]
            else:
                triangles_reduced = None
            triangles_reduced = comm.bcast(triangles_reduced, root=0)
            triangles_reduced = np.array(triangles_reduced)
            number_of_triangles = len(triangles)
            
            del triangles
            del up

            loc_tri = distribute(triangles_reduced, rank, size)
            
            dos_loc, nd_loc = bzi.linear_triangulation(nbands, loc_tri, len(E), E,
                                                       Freq, func, rank) 
            
            DOS = np.zeros(len(dos_loc),dtype=dos_loc.dtype)
            ND = np.zeros(len(nd_loc),dtype=nd_loc.dtype)
            comm.Allreduce(dos_loc,DOS,op=MPI.SUM)
            comm.Allreduce(nd_loc,ND,op=MPI.SUM)
            
            DOS *= (1/number_of_triangles/3)/moire_area 
            ND *= (1/number_of_triangles/3)/moire_area *3*self.natom/nbands
            
            if rank == 0:
                f = open(output_file,'w+')
                for i in range(len(E)):
                    f.writelines("%f\t\t%f\t%f\n"%(E[i],DOS[i],ND[i]))
                f.close()
                print("\nWritten DOS and Number density to %s"%(output_file), flush=True)
        
        elif (method=='gaussian') or (method=='lorentzian'): 
            if method=='gaussian':
                mm = 1    
            elif method=='lorentzian':
                mm = 2
            weights = np.loadtxt(q_grid_file,usecols=(3,))
            data = h5py.File(data_file,'r')
            global_keys = list(data.keys())
            local_keys = distribute(global_keys,rank,size)
            weights_local = np.empty(len(local_keys),dtype='f')
            for i,keys in enumerate(local_keys):
                weights_local[i] = weights[int(keys)-1]
            dos_loc = np.zeros(np.shape(E))
            for group in local_keys:
                Freq = data[group]['eigenvalues'][:]
                dos_loc += bzi.standard_integration(method=mm,weight=weights_local,
                                                    freq=Freq,width=width,energy=E,
                                                    nenergy=len(E),natom=self.natom,
                                                    nbands=nbands)
            DOS = np.zeros(len(dos_loc),dtype=dos_loc.dtype)
            comm.Allreduce(dos_loc,DOS,op=MPI.SUM)
            if rank==0:
                f = open(output_file,'w+')
                for i in range(len(E)):
                    f.writelines("%f\t\t%f\n"%(E[i],DOS[i]))
                f.close()
                print("\nWritten DOS to %s"%(output_file), flush=True)
        else:
            if rank==0:
                print("Method not recognised",flush=True)
                print("Recognised methods are: ",flush=True)
                print(" 1. linear_triangulation",flush=True)
                print(" 2. gaussian",flush=True)
                print(" 3. lorentzian",flush=True)
                print(" ")
            exit()
            
        return    



    def compute_polarization(self, evec):
        """
            Computes the polarization of a particular mode
            
            Input
            -----
                evec: Eigenvector (numpy array of length 3*natom)
            Output
            ------

                sx,sy,sz: polarization in x,y and z

            Theory
            ------
                s_x = 2*[Re(e_y)*Im(e_z) - Im(e_y)*Re(e_z)]
                s_y = 2*[Re(e_z)*Im(e_x) - Im(e_z)*Re(e_x)]
                s_z = 2*[Re(e_x)*Im(e_y) - Im(e_x)*Re(e_y)]
        """ 
        evec = evec.reshape(-1,3)
        ex_r,ex_i = evec[:,0].real, evec[:,0].imag
        ey_r,ey_i = evec[:,1].real, evec[:,1].imag
        ez_r,ez_i = evec[:,2].real, evec[:,2].imag
        s_x = 2*(ey_r@ez_i - ey_i@ez_r)
        s_y = 2*(ez_r@ex_i - ez_i@ex_r)
        s_z = 2*(ex_r@ey_i - ex_i@ey_r)
        return s_x, s_y, s_z


    def chirality(self,
                  ph_data_file, 
                  output_file_name,
                  nbands):
        """
            Writes the chirality data for all q-points

            ph_data_file :: phonon data output file containing the eigenvectors

        """

        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()


        f = h5py.File(ph_data_file,'r')
        groups = list(f.keys())

        nqpt = len(groups)

        loc_groups = distribute(groups,rank,size)

        out_f = h5py.File(output_file_name,'w',driver='mpio',comm=comm)
        dset_out = out_f.create_dataset('chirality',((nqpt,nbands,3)),dtype=np.float64)
        
        for lg in loc_groups:
            for nu in range(nbands):
                evec_q = f[lg]['evec'][nu]
                sx,sy,sz = self.compute_polarization(evec_q)
                dset_out[int(lg)-1,nu,:] = np.array([sx,sy,sz],dtype=np.float64)
            print("Rank %d finished computations at q-point %d"%(rank,int(lg)),flush=True)
        f.close()
        out_f.close()

        return
