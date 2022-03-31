from pymelecutil import plot_figures

label=(r'$\Gamma $','M', 'K', r'$\Gamma $')
en_range = [-.250,.250]
plt = plot_figures(dpi=500,en_range=en_range)
plt.plot_band_structure(data_file='/home/mshinjan/Electronic_bands_data/tblg_0.95/bands_0.95_E0.hdf5',label=label,
                        nbands=14284,nkpt=70,kfile='/home/mshinjan/Electronic_bands_data/tblg_0.95/kpt_bands.dat',
                        save=False)
#plt.plot_band_structure(data_file='../Examples/tblg_21.8/bands_21.8_E0.hdf5',label=label, 
#                        nbands=28, nkpt=300, kfile='../Examples/tblg_21.8/k_points.dat',
#                        en_range=en_range,save=False)
