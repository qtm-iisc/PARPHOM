import matplotlib.pyplot as plt
import numpy as np
import h5py
from pyphutil.pyphutil import moire_phonon_utils
from mpl_toolkits.axes_grid1 import make_axes_locatable


mph = moire_phonon_utils(lammps_input_file='/home/mshinjan/Phonons/MoirePhonons/Force_Constant_Generation/Examples/tbwse2_21.8/in.phonon')

mph.read_lammps()

print(mph.rec_lat)


chiral_file = 'chirality_dos_21.8_tbwse2.hdf5'
dos_file = 'phbands_dos_21.8.hdf5'



f_dos = h5py.File(dos_file,'r')
f_chi = h5py.File(chiral_file,'r')

dos_keys = list(f_dos.keys())

q_vec = []

for i,q in enumerate(dos_keys):
    q_vec.append(np.matmul(mph.rec_lat.T, f_dos[q]['q_vec'][:]))

q_vec = np.array(q_vec)
chiral = f_chi['chirality'][:]

from mpl_toolkits.axes_grid1 import make_axes_locatable

for i in range(chiral.shape[1]):
    fig,axs = plt.subplots(1,3,figsize=(27*0.75,9*0.75))#,constrained_layout=True)
    im0 = axs[0].scatter(q_vec[:,0],q_vec[:,1],c=chiral[:,i,0],vmax=1.0, vmin=-1.0,cmap='coolwarm')
    im1 = axs[1].scatter(q_vec[:,0],q_vec[:,1],c=chiral[:,i,1],vmax=1.0, vmin=-1.0,cmap='coolwarm')
    im2 = axs[2].scatter(q_vec[:,0],q_vec[:,1],c=chiral[:,i,2],vmax=1.0, vmin=-1.0,cmap='coolwarm')
    axs[0].set_aspect('equal')
    axs[1].set_aspect('equal')
    axs[2].set_aspect('equal')
    divider0 = make_axes_locatable(axs[0])
    cax0 = divider0.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im0, cax=cax0)
    divider1 = make_axes_locatable(axs[1])
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im1, cax=cax1)
    divider2 = make_axes_locatable(axs[2])
    cax2 = divider2.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im2, cax=cax2)
    plt.tight_layout()
    plt.savefig("Chiral_%03d.png"%i)
    plt.close(fig)
    print("%3d done"%i,flush=True)
