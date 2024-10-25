"""
    Written by Shinjan Mandal as an utility for PARPHOM
    Github repo: https://github.com/qtm-iisc/PARPHOM
    Email : mshinjan@iisc.ac.in
"""
import matplotlib.pyplot as plt
import numpy as np
import h5py

class plot():
    """
        Class for plotting figures 
    """

    def __init__(self, 
                 en_range=None,
                 dpi=None):
    
        self.en_range = en_range
        self.dpi = dpi

    def band_structure(self,
                       data_file,
                       nbands,
                       nqpt,
                       qfile,
                       label,
                       save=True,
                       savefile='bandstructure.png',
                       closed_loop=True):
        """
            Plot the band structure of a generated phonon spectrum
        """
        
        with open(qfile,'r') as f:
            for nodes in f:
                pass
        nodes = nodes.split()
        nodes = np.array([int(i) for i in nodes])
        index_ = [i for i in range(nqpt)]
        f = h5py.File(data_file,'r')
        
        groups = list(f.keys())

        if closed_loop:
            eigvals = np.empty((nqpt+1, nbands))
        else:
            eigvals = np.empty((nqpt,nbands))

        counter = 0
        for group in groups:
            eigvals[counter] = f[group]['eigenvalues'][:]
            counter += 1
        if (closed_loop):
            eigvals[nqpt] = eigvals[0]
            counter += 1
        x = [i for i in range(counter)]
        for i in range(nbands):
            plt.plot(x,eigvals[:,i],c='b')
        plt.xticks(nodes,label,fontsize=22)
        plt.ylabel(r'Energy (cm$^{-1}$)', size=22)
        plt.yticks(fontsize=20)
        if self.en_range!=None:
            plt.ylim(self.en_range[0],self.en_range[1])
        for n in range(len(nodes)):
            plt.axvline(x=nodes[n],linewidth=1,color='k')
        plt.axhline(y=0,linewidth=1,color='k')
        plt.xlim(nodes[0],nodes[len(nodes)-1])
        if save:
            if self.dpi!=None:
                plt.savefig(savefile,bbox_inches='tight',dpi=self.dpi)
            else:
                plt.savefig(savefile,bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()
        plt.clf()
        return 




    def dos(self,
            dos_file,
            save=True,
            savefile='dos.png',
            max_dos=None):
        """
            Plot the DOS
        """

        dos_data = np.loadtxt(dos_file)
        plt.fill(dos_data[:,0],dos_data[:,1],color='slategrey')
        plt.xlabel(r'Energy (cm$^{-1}$)', size=22)
        plt.ylabel(r'DOS(cm$^{-1}$nm$^{-2}$)',fontsize=22)
        plt.tick_params(axis='both', which='major', labelsize=20)
        if max_dos!=None:
            plt.ylim(0,max_dos)
        else:
            plt.ylim(0,)
        if self.en_range!=None:
            plt.xlim(self.en_range[0],self.en_range[1])
        if save:
            if self.dpi!=None:
                plt.savefig(savefile,bbox_inches='tight',dpi=self.dpi)
            else:
                plt.savefig(savefile,bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()
        plt.clf()
        return

