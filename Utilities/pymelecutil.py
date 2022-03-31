"""
   Written by Shinjan Mandal as an utlity for Moire Electronic Bands
   Github repo: https://github.com/ShinjanM/MoireElectronicBands
   Email : mshinjan@iisc.ac.in   

"""

import numpy as np
import spglib
import matplotlib.pyplot as plt
from mpi4py import MPI
import h5py

class moire_electron_utils():
    """
        Utility routines for Moire Electronic Bands
    """
    def __init__(self,
                 file_name=None,
                 file_location=None,
                 natom=None,
                 at_types=None,
                 at_style=None):
        self.file_name = file_name
        self.file_location = file_location
        self.natom = natom
        self.at_types = at_types
        self.at_style = at_style
        self.at_id = None
        self.mol_id = None
        self.mass = None
        self.lat = None
        self.rec_lat = None
        self.real_pos = None
        self.crys_pos = None
        self.neigh_list = None

    def read_lammps(self):
        """
            Reads the LAMMPS output file 
            Inputs needed for the lammps_data class:
                file_name - Name of the LAMMPS file
                file_location - Location of the LAMMPS file
                natom - Number of atoms in the system
                at_types - Number of atom types
                at_style - Atomic style used in the LAMMPS data file

        """
        if self.file_location[len(self.file_location)-1]!='/':
            self.file_location += '/'
        fname = self.file_location+self.file_name 
        f = open(fname,'r')
        lines = f.readlines()
        f.close()
        for i in range(len(lines)):
            if "xlo xhi" in lines[i]:
                lx = eval(lines[i].split()[1]) - eval(lines[i].split()[0])
            if "ylo yhi" in lines[i]:
                ly = eval(lines[i].split()[1]) - eval(lines[i].split()[0])
            if "zlo zhi" in lines[i]:
                lz = eval(lines[i].split()[1]) - eval(lines[i].split()[0])
            if "xy xz yz" in lines[i]:
                xy = eval(lines[i].split()[0])
                xz = eval(lines[i].split()[1])
                yz = eval(lines[i].split()[2])
            if "Atoms" in lines[i]:
                if (self.at_style == 'atomic'):
                    self.at_id =  np.array([eval(lines[j].split()[0]) for j in 
                                                range(i+2, i+self.natom+2)])
                    self.mol_id = np.array([eval(lines[j].split()[1]) for j in 
                                            range(i+2, i+self.natom+2)])
                    self.real_pos = np.array([[eval(lines[j].split()[2]),
                                               eval(lines[j].split()[3]),
                                               eval(lines[j].split()[4])]
                                               for j in range(i+2, i+self.natom+2)])
                elif (self.at_style == 'molecular'):
                    self.at_id =  np.array([eval(lines[j].split()[0]) for j in 
                                                range(i+2, i+self.natom+2)])
                    self.mol_id = np.array([eval(lines[j].split()[2]) for j in 
                                                range(i+2, i+self.natom+2)])
                    self.real_pos = np.array([[eval(lines[j].split()[3]),
                                               eval(lines[j].split()[4]),
                                               eval(lines[j].split()[5])]
                                               for j in range(i+2, i+self.natom+2)])
            if "Masses" in lines[i]:
                self.mass = np.array([eval(lines[j].split()[1]) for j in 
                                    range(i+2,i+self.at_types+2)])
        
        
        """

            For determining the lattice the convention used is:
            
                                          | ax, bx, cx |
            lattice vector is stored as : | ay, by, cy |
                                          | az, bz, cz |

            where the lattice vectors are: a, b and c

                                          | Ga1,  Gb1,  Gc1 |
            reciprocal lattice vector is: | Ga2,  Gb2,  Gc2 |
                                          | Ga3,  Gb3,  Gc3 |

        """
        
        self.lat = np.array([[ lx, xy, xz],
                             [  0, ly, yz],
                             [  0,  0, lz]])
        self.rec_lat = (np.linalg.inv(self.lat)).transpose()*2*np.pi
       
        self.crys_pos = np.array([np.linalg.solve(self.lat,self.real_pos[i])
                                        for i in range(self.natom)])
        
        return


    def get_kpt_bz(self, 
                   mesh, 
                   use_symmetry=True, 
                   GAMMA=True,
                   output_file_name='kpt_dos.dat',
                   write_full_grid=True,
                   full_grid_file_name='kpt_full.dat'):
        """
            Generates the points on the k-grid mesh you want to sample. 
    
            Input:  Enter the mesh as a numpy array
                    Example: mesh = np.array([5,5,1]), which denotes 
                             a 5X5X1 mesh
                    
                    GAMMA: True for gamma centred mesh, False otherwise  
            
            Output: 
                    k_points_crys: numpy array of the k points to be sampled 
                                   in crystal coordinates 
                    weights: the weights corresponding to each k-point
        """
        
        cell = (self.lat, self.crys_pos, self.mol_id)
        if GAMMA==True:
            # Γ centered Mesh
            mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0,0,0])
            if use_symmetry:
                u, inverse = np.unique(mapping, return_inverse=True)
                k_points_crys = (grid[np.unique(u)])/np.array(mesh,dtype=float) 
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
                k_points_crys = grid/mesh
        else:
            # Γ shifted Mesh
            mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0.5,0.5,0])
            if use_symmetry:
                u, inverse = np.unique(mapping, return_inverse=True)
                k_points_crys = (grid[np.unique(u)]+[0.5,0.5,0.0])/np.array(mesh,dtype=float)
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
                k_points_crys = (grid+np.array([0.5,0.5,0.0]))/mesh
        nkpt = mesh[0]*mesh[1]*mesh[2]
        if use_symmetry:
            weights = np.array([0 for i in range(len(k_points_crys))])
            for i in inverse:
                weights[i] += 1
        else:
            weights = np.array([1 for i in range(len(k_points_crys))])
        weights = weights/nkpt
        f = open(output_file_name,'w')
        for i in range(len(k_points_crys)):
            f.write('%+.6f\t%+.6f\t%+.6f\t\t%+.8f\n'%(k_points_crys[i,0],
                                             k_points_crys[i,1],
                                             k_points_crys[i,2],
                                             weights[i]))
        return



    def get_kpt_bands(self, points,output_file='k_points.dat'):
        """
            Generates a list of points along the path specified for plotting the 
            band structure

            Input: 
                  points: list of k-points along which the bands are plotted
            
                  The input style of points is - 
                  [..., k(i), no of points between k(i) 
                                        and k(i+1), k(i+1), ...] 
                  where k(n) is the n-th k-point, in fractional coordinates
                                            
        
            Output: 
                  k_pt: list of the k-point path in fractional coordinates
        """
        k_pt, nodes = [], [0,]
        for i in range(0,len(points)-1,2):
            ki = points[i]
            kiplus1 = points[i+2]
            num = points[i+1]
            nodes.append(num+nodes[len(nodes)-1])
            for i in range(num):
                k_pt.append(ki+((kiplus1-ki)*i/num))
        if (np.allclose(points[0],points[len(points)-1])==False):
            k_pt.append(points[len(points)-1])
        k_pt = np.array(k_pt)
        f = open(output_file,'w')
        for i in range(len(k_pt)):
            f.write('%+.6f\t%+.6f\t%+.6f\n'%(k_pt[i,0],k_pt[i,1],k_pt[i,2]))
        for i in range(len(nodes)):
            f.write('%d\t'%nodes[i])
        return


    def read_full_grid(self, mesh, full_grid_file):
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
                          method='linear triangulation',
                          k_grid_file = None,
                          width = None,
                          full_k_grid_file=None,
                          output_file='dos.dat'):
        """
            Computes the density of states of a system
        """

        from bz_integration import bz_integration as bzi
        from distribute_lists import distribute

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()  
        size = comm.Get_size()
        

        if method=='linear triangulation':

            from scipy.spatial import Delaunay

            moire_area = abs(np.linalg.det(self.lat[0:2,0:2]))*0.01 # in nm^2
            BZ_area = abs(np.linalg.det(self.rec_lat[0:2,0:2]))*100 # in nm^-2

            grid, mapping = self.read_full_grid(mesh, full_k_grid_file)
            
            Freq = np.empty((len(np.unique(mapping)),nbands))
            data = h5py.File(data_file,'r')
            counter = 0
            for group in data.keys():
                ds_data = data[group]['eigenvalues']
                Freq[counter] = ds_data[:]
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
                print("Rank %d appended %d/%d triangles."%(rank,i+1,len(local_triangles)),flush=True)
            comm.Barrier()

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
            E = np.arange(en_range[0],en_range[1],en_spacing)
            func = np.ones(np.shape(Freq),dtype=np.float64)
            
            dos_loc, nd_loc = bzi.linear_triangulation(nbands, loc_tri, len(E), E,
                                                       Freq, func, rank) 
            
            Master_list_global_dos = []
            Master_list_global_dos = comm.gather(dos_loc, root=0)
            Master_list_global_ndos = []
            Master_list_global_ndos = comm.gather(nd_loc, root=0)

            if rank == 0:
                DOS = np.zeros(len(dos_loc))
                ND = np.zeros(len(nd_loc))
                for i in range(size):
                    DOS += Master_list_global_dos[i]
                    ND += Master_list_global_ndos[i]
                DOS *= (2/number_of_triangles/3)/moire_area 
                ND *= (2/number_of_triangles/3)/moire_area *self.natom/nbands

                f = open(output_file,'w+')
                for i in range(len(E)):
                    f.writelines("%f\t\t%f\t%f\n"%(E[i],DOS[i],ND[i]))
                f.close()
                print("\nWritten DOS to %s"%(output_file), flush=True)

        return        
        

    def interlayer_seperation(self,
                              no_of_supercells,
                              box_dim,
                              image_size=5000,
                              number_of_ticks = 5,
                              vmin = 3.35,
                              vmax = 3.65,
                              dpi = 300,
                              output_file='interlayer_seperation.pdf',
                              transparent=True):
        """
            Computes the interlayer seperation for a given moire structure.
            
            Input:
            -----
                no_of_supercells (int) == Number of supercells in each direction 
                                          to be scanned
  
                box_dim (numpy array of dim(2,2)) == Dimensions of the region to 
                                                     scan.
        
                                                    | x_min,   x_max |
                        The syntax of box_dim is:   |                |
                                                    | y_min,   y_max |
                
                image_size (int) == Number of real space points in the image
                (default = 5000)

                number_of_ticks (int) == Number of ticks in each axis.
                 (default = 5)           The ticks are equally spaced among 
                                         axis min values and axis max value.

                vmin (float) == Minumum value of height to be shown in the cmap
                vmax (float) == Maxiumum value of height to be shown in the cmap
                
                dpi (int) == dpi for matplotlib figure
            
            Output:
            -------

                output_file : Name of the output file where the 
                              interlayer seperation data is plotted.
        
        """

        from scipy.interpolate import griddata
        from matplotlib import cm

        X1,Y1,Z1 = [],[],[]
        X2,Y2,Z2 = [],[],[]

        crys = np.split(self.crys_pos, 2)
        
        for item in crys[0]:
            for i in range(-no_of_supercells, no_of_supercells+1):
                for j in range(-no_of_supercells, no_of_supercells+1):
                    pt_crys = np.array([item[0]+i, item[1]+j, item[2]])
                    pt_real = np.matmul(self.lat.T, pt_crys)
                    if (pt_real[0]>=box_dim[0,0] and pt_real[0]<=box_dim[0,1]) and \
                       (pt_real[1]>=box_dim[1,0] and pt_real[1]<=box_dim[1,1]):
                           X1.append(pt_real[0])
                           Y1.append(pt_real[1])
                           Z1.append(pt_real[2])

        for item in crys[1]:
            for i in range(-no_of_supercells, no_of_supercells+1):
                for j in range(-no_of_supercells, no_of_supercells+1):
                    pt_crys = np.array([item[0]+i, item[1]+j, item[2]])
                    pt_real = np.matmul(self.lat.T, pt_crys)
                    if (pt_real[0]>=box_dim[0,0] and pt_real[0]<=box_dim[0,1]) and \
                       (pt_real[1]>=box_dim[1,0] and pt_real[1]<=box_dim[1,1]):
                           X2.append(pt_real[0])
                           Y2.append(pt_real[1])
                           Z2.append(pt_real[2])

        box_x = np.linspace(box_dim[0,0]+1.5, box_dim[0,1]-1.5, image_size) 
        box_y = np.linspace(box_dim[1,0]+1.5, box_dim[1,1]-1.5, image_size)

        X,Y = np.meshgrid(box_x, box_y)

        box_x = np.linspace(0,image_size,number_of_ticks)
        box_y = np.linspace(0,image_size,number_of_ticks)

        L1 = griddata((X1,Y1),Z1,(X,Y),method='cubic')
        L2 = griddata((X2,Y2),Z2,(X,Y),method='cubic')

        ils = L2-L1

        cmap = plt.cm.get_cmap('inferno')
            
        plt.imshow(ils,cmap=cmap, vmin=vmin, vmax=vmax)

        cbar = cm.ScalarMappable(cmap=cmap)
        if vmin!=None and vmax!=None:
            cbar.set_clim(vmin=vmin, vmax=vmax)

        xticks = np.linspace(box_dim[0,0]+1.5,box_dim[0,1]-1.5,number_of_ticks)
        for i in range(len(xticks)):
            xticks[i] = round(xticks[i],2)
        yticks = np.linspace(box_dim[1,0]-1.5,box_dim[1,1]+1.5,number_of_ticks)
        for i in range(len(yticks)):
            yticks[i] = round(yticks[i],2)

        plt.xticks(box_x,xticks,fontsize=14)
        plt.yticks(box_y,yticks,fontsize=14)

        plt.xlabel(r'X ($\AA$)',size=15)
        plt.ylabel(r'Y ($\AA$)',size=15)
        Cbar = plt.colorbar(cbar)
        Cbar.ax.tick_params(labelsize=12)
        Cbar.set_label(r'Interlayer distance ($\AA$)',size=15,rotation=270,labelpad=40)
        
        plt.tight_layout()
        plt.savefig(output_file, format='pdf', dpi=dpi, transparent=transparent)
        plt.clf()

        return



    def strain(self, 
               equib_latcon, 
               neighbors,
               heading=None, 
               clim=None, 
               save=False, 
               dpi=300,
               plot_limits=None,
               output_file_name='Strain.jpg'):
        from neighbor import neighbor_list as n
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from mpl_toolkits.axes_grid1 import ImageGrid

        if self.neigh_list==None:
            self.neigh_list = n.create_neigh_list(self.real_pos, self.lat)
        strain = n.calculate_strain(self.neigh_list, equib_latcon)
        #fig,ax = plt.subplots(1,2)
        fig = plt.figure(figsize=(18,6))
        grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,2),
                 axes_pad=0.25,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )
        ax = []
        for j in grid:
            ax.append(j)
            j.set_aspect('equal')
            j.set_xlabel(r"x ($\AA$)", fontsize = 18)
            j.set_ylabel(r"y ($\AA$)", fontsize = 18)

        lines = []
        c = []
        for i in range(self.natom//2):
            for ix in range(-neighbors,neighbors+1):
                for iy in range(-neighbors, neighbors+1):
                    trans = np.matmul(self.lat[0:2,0:2], np.array([ix,iy]))
                    p0 = self.real_pos[i,:2] +  trans
                    p1 = self.neigh_list[i,0,4:6] +  trans
                    p2 = self.neigh_list[i,1,4:6] +  trans
                    p3 = self.neigh_list[i,2,4:6] +  trans
                    lines.extend([[p0,p1],[p0,p2],[p0,p3]])
                    c.extend([strain[i,0], strain[i,1], strain[i,2]])
        line_segments = LineCollection(np.array(lines), linestyles='solid',
                                       cmap='coolwarm',lw=0.5)
        c = np.array(c)
        line_segments.set_array(c)
        ax[0].add_collection(line_segments)
        axcb0 = ax[0].cax.colorbar(line_segments)
        axcb0.ax.tick_params(labelsize=15)
        if clim!=None:
            axcb0.mappable.set_clim(clim[0],clim[1])
        ax[0].autoscale()

        lines = []
        c = []
        for i in range(self.natom//2,self.natom):
            for ix in range(-neighbors,neighbors+1):
                for iy in range(-neighbors, neighbors+1):
                    trans = np.matmul(self.lat[0:2,0:2], np.array([ix,iy]))
                    p0 = self.real_pos[i,:2] +  trans
                    p1 = self.neigh_list[i,0,4:6] +  trans
                    p2 = self.neigh_list[i,1,4:6] +  trans
                    p3 = self.neigh_list[i,2,4:6] +  trans
                    lines.extend([[p0,p1],[p0,p2],[p0,p3]])
                    c.extend([strain[i,0], strain[i,1], strain[i,2]])
        line_segments = LineCollection(np.array(lines), linestyles='solid',
                                       cmap='coolwarm',lw=0.5)
        c = np.array(c)
        line_segments.set_array(c)
        ax[1].add_collection(line_segments)
        axcb1 = ax[1].cax.colorbar(line_segments)
        axcb1.ax.tick_params(labelsize=15)
        if clim!=None:
            axcb1.mappable.set_clim(clim[0],clim[1])
        ax[1].autoscale()


        a1 = self.lat[:,0][0:2]
        a2 = self.lat[:,1][0:2]
        orig = np.array([0,0])

        pc_boundary = LineCollection([[orig,a1],[orig,a2],[a1,a1+a2],[a2,a1+a2]],
                                     color='gray',lw=4.5)
        ax[0].add_collection(pc_boundary)
        pc_boundary = LineCollection([[orig,a1],[orig,a2],[a1,a1+a2],[a2,a1+a2]],
                                     color='gray',lw=4.5)
        ax[1].add_collection(pc_boundary)

        for j in grid:
            if plot_limits==None:
                j.set_xlim(-self.lat[0,0]*neighbors/2+self.lat[0,0]*0.75,
                            self.lat[0,0]*neighbors/2+self.lat[0,0]*0.75)
                j.set_ylim(-self.lat[1,1]*neighbors/2+self.lat[1,1]*0.75,
                            self.lat[1,1]*neighbors/2+self.lat[1,1]*0.75)
            else:
                j.set_xlim(-plot_limits[0],plot_limits[0])
                j.set_ylim(-plot_limits[1],plot_limits[1])
            j.tick_params(axis='both',labelsize=15)
            p1 = 0*a1+0*a2
            p2 = 1/3*a1+1/3*a2
            p3 = 2/3*a1+2/3*a2
            p4 = a1+a2
            j.plot(p1[0],p1[1],'o',c='r')
            j.plot(p2[0],p2[1],'o',c='g')
            j.plot(p3[0],p3[1],'o',c='k')
            j.plot(p4[0],p4[1],'o',c='r')

        

        if heading!=None:
            fig.suptitle(heading, size=22)
        if save:
            plt.savefig(output_file_name, dpi=dpi, bbox_inches='tight')
            plt.clf()
        else:
            plt.show()
        return












class plot_figures():
    """
        Plot the band structures, DOS, velocity and other relevant data extracted
    """
    

    def __init__(self,
                 fig_size = None,
                 dpi = None,
                 en_range=None):
        self.fig_size = fig_size
        self.dpi = dpi
        self.en_range = en_range
        return

    def plot_band_structure(self, data_file, nbands, nkpt, kfile, label, 
                            save=True,savefile='bandstruct.png', closed_loop=True):
        
        with open(kfile,'r') as f:
            for nodes in f:
                pass
        nodes = nodes.split()
        nodes = np.array([int(i) for i in nodes])
        index_ = [i for i in range(nkpt)]
        data = h5py.File(data_file,'r')
        if closed_loop:
            eigvals = np.empty((nkpt+1, nbands))
        else:
            eigvals = np.empty((nkpt,nbands))

        counter = 0
        for group in data.keys():
            ds_data = data[group]['eigenvalues']
            eigvals[counter] = ds_data[:]
            counter += 1
        if (closed_loop):
            eigvals[nkpt] = eigvals[0]
            counter += 1
        x = [i for i in range(counter)]
        for i in range(nbands):
            plt.plot(x,eigvals[:,i],c='b')
        plt.xticks(nodes,label,fontsize=18)
        plt.ylabel('Energy (eV)', size=15)
        if self.en_range!=None:
            plt.ylim(self.en_range[0],self.en_range[1])
        for n in range(len(nodes)):
            plt.axvline(x=nodes[n],linewidth=1,color='k')
        plt.axhline(y=0,linewidth=1,color='k')
        plt.xlim(nodes[0],nodes[len(nodes)-1])
        if save:
            plt.savefig(savefile,dpi=self.dpi)
        else:
            plt.show()
        plt.clf()
        return
