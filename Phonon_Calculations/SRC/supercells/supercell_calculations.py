import h5py
import numpy as np


class supercell_calculations():

    def __init__(self, lammps_input, fc_file, q_file):
        
        self.lammps_input = lammps_input
        self.fc_file = fc_file
        self.q_file = q_file
    
    def read_force_constant_file(self):
        f = h5py.File(fc_file,'r')
        fc = f['force_constants'][:]
        p2s_map = f['p2s_map'][:]
        p_indices = f['p_indices'][:]
        self.force_constants = fc
        self.prim2sup_map = p2s_map
        self.prim_inds = p_indices
        return 

    def read_q(self):
        q_pts = np.genfromtext(self.q_file, skip_footer=1)
        return q_pts

    
    def read_lammps(self, show_log=False):
        """
            Read lammps data file and extract data
        """
        from lammps import lammps
        cmd_list = ['-log', 'none']
        if not show_log:
            cmd_list += ['-echo', 'none', '-screen', 'none']
        lmp = lammps(cmdargs=cmd_list)
        lmp.file(self.lammps_input)
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


    def create_dynamical_matrix(self, q):
        """
            Create the dynamical matrix for a supercell
        """
        return
