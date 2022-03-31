from pymelecutil import plot_figures

label=(r'$\Gamma $','M', 'K', r'$\Gamma $', r'K$^{p}$')
en_range = None #[-.250,.250]
plt = plot_figures(dpi=500,en_range=en_range)
plt.plot_band_structure(data_file='../Examples/tbmos2_21.8/phbands_21.8_mos2_KKp.hdf5',
                        label=label, 
                        nbands=126, 
                        nkpt=181, 
                        kfile='../Examples/tbmos2_21.8/qpt_bands_KKp.dat',
                        save=False,
                        with_chirality=True,
                        closed_loop=False,
                        chiral_file='/home/mshinjan/Phonons/Phonon_ang_momentum/test_file_bands')
