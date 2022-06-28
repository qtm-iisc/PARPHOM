import h5py
import numpy as np
from lammps import lammps
import matplotlib.pyplot as plt

amp = 1e-4

cmd_list = ['-log', 'none','-echo', 'none', '-screen', 'none']

lmp = lammps(cmdargs=cmd_list)

lmp.file('../../../Force_Constant_Generation/Examples/tblg_21.8/in.phonon')
lmp.command('run 0')
natom = lmp.get_natoms()
id_ = lmp.extract_atom("id",0)
id_ = np.array([id_[i] for i in range(natom)],dtype=int)
xp = lmp.extract_atom("x",3)
ref = np.array([[xp[i][0], xp[i][1], xp[i][2]] for i in range(natom)], dtype=float)
e_ref = lmp.get_thermo('etotal')
print("Total energy = ", e_ref)


ph_vec_file_phonopy = '../Phonopy_Comparision/Phonopy_Comparision_21_8.h5'
Phonopy_units = 4.13567*8.065610

f = h5py.File(ph_vec_file_phonopy, 'r')
groups = list(f.keys())
groups = groups[:1]
disp_freq_phonopy = np.zeros((len(groups),3*natom),dtype=float)
eigvals_phonopy = np.zeros((len(groups),3*natom),dtype=float)
ph_vec_gamma_phonopy = f[groups[0]]['evec'][:]
eigvals_phonopy[0] = f[groups[0]]['eigenvalues'][:]*Phonopy_units
for i in range(len(ph_vec_gamma_phonopy)):
    for n in range(natom):
        lmp.command('set atom {} x {} y {} z {}'.format(id_[n], 
                                                ref[n,0]+amp*ph_vec_gamma_phonopy[i,3*n+0].real,
                                                ref[n,1]+amp*ph_vec_gamma_phonopy[i,3*n+1].real,
                                                ref[n,2]+amp*ph_vec_gamma_phonopy[i,3*n+2].real))
    lmp.command('run 0')
    e_mode = lmp.get_thermo('etotal')
    freq = np.sqrt(abs(((e_mode-e_ref)*1000*2*4.13567)/(amp**2 * 12.001)))*8.065610
    print(freq, eigvals_phonopy[0,i],flush=True)
    disp_freq_phonopy[0,i] = (eigvals_phonopy[0,i]-freq)/freq * 100

ph_file = 'phbands_21.8.hdf5'
f = h5py.File(ph_file, 'r')
groups = list(f.keys())
groups = groups[:1]
disp_freq = np.zeros((len(groups),3*natom),dtype=float)
eigvals = np.zeros((len(groups),3*natom),dtype=float)
ph_vec_gamma = f[groups[0]]['evec'][:]
eigvals[0] = f[groups[0]]['eigenvalues'][:]
for i in range(len(ph_vec_gamma)):
    for n in range(natom):
        lmp.command('set atom {} x {} y {} z {}'.format(id_[n], 
                                                ref[n,0]+amp*ph_vec_gamma[i,3*n+0].real,
                                                ref[n,1]+amp*ph_vec_gamma[i,3*n+1].real,
                                                ref[n,2]+amp*ph_vec_gamma[i,3*n+2].real))
    lmp.command('run 0')
    e_mode = lmp.get_thermo('etotal')
    freq = np.sqrt(abs(((e_mode-e_ref)*1000*2*4.13567)/(amp**2 * 12.001)))*8.065610
    print(freq, eigvals[0,i],flush=True)
    disp_freq[0,i] = (eigvals[0,i]-freq)/freq * 100


x = np.array([i for i in range(1,3*natom+1)])

fig,axs = plt.subplots(figsize=(16*0.9,9*0.9),constrained_layout=True)

axs.bar(x,height=disp_freq[0],width=0.6,label=r"$\Delta_{\mathbf{q}=\Gamma,\nu}$ current work",zorder=2,color='b',alpha=0.6)
axs.bar(x,height=disp_freq_phonopy[0],width=0.9,alpha=0.8,zorder=1,label=r"$\Delta_{\mathbf{q}=\Gamma,\nu}$ Phonopy",color='g')
#axs.bar(x,height=disp_freq_phonopy[0]-disp_freq[0],width=0.45,label=r"Difference in error",zorder=3,color='r',alpha=1.0)
axs.set_xlabel(r"$\nu$ (Mode Number)",fontsize=22)
axs.set_ylabel(r"\% Error",fontsize=22)
axs.tick_params(axis='both',which='major',labelsize=18)
axs.legend(loc='best',fontsize=16,fancybox=True,shadow=True)
axs.set_xlim(5.5,85)
axs.set_ylim(-0.5,6)
plt.show()
