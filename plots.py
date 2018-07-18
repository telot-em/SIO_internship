import func_plot
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


plot_path=None #'/home/ebent/plots/figs/'

f = Dataset('aviso_save.nc', 'r', format='NETCDF4')
mean_aviso = f.variables['mean_aviso'][:]
std_aviso = f.variables['std_aviso'][:]
Ya = f.variables['Ya'][:]
Xa = f.variables['Xa'][:]


f2 = Dataset('model_save.nc', 'r', format='NETCDF4')
mean_model = f2.variables['mean_model'][:]
std_model = f2.variables['std_model'][:]

f3 = Dataset('diff_save.nc', 'r', format='NETCDF4')
diff = f3.variables['diff_weight'][:]
mean_diff = f3.variables['mean_diff_weight'][:]


# AVISO DATA
func_plot.polar_map(np.linspace(0, 0.52, 50), Xa, Ya, std_aviso, 'AVISO data, std', 'Standard deviation [m]', plt.cm.jet, plot_path, title_save='aviso_std.png')

func_plot.polar_map(np.linspace(-2, 2, 50), Xa, Ya, mean_aviso, 'AVISO data, mean', 'MADT [m]', plt.cm.jet, plot_path, 'aviso_mean.png')


# MODEL OUTPUT
func_plot.polar_map(np.linspace(0, 0.52, 50), Xa, Ya, std_model, 'Model output, std', 'Standard deviation [m]', plt.cm.jet, plot_path, 'model_std.png')

func_plot.polar_map(np.linspace(0, 0.52, 50), Xa, Ya, mean_model, 'Model output, mean', 'MADT [m]', plt.cm.jet, plot_path, 'model_mean.png')




# Difference between data and model
func_plot.polar_map(np.linspace(-1, 1, 50), Xa, Ya, diff, 'Difference of mean SSH between data and model', 'MADT [m]', plt.cm.bwr, plot_path, 'diff.png')




# Add the offset to model
func_plot.polar_map(np.linspace(-2,2,50), Xa, Ya, mean_model+mean_diff, 'Mean SSH of model + offset', 'MADT [m]', plt.cm.jet, plot_path, 'offset+model.png')

# Diff - offset
func_plot.polar_map(np.linspace(-0.5,0.5,50), Xa, Ya, diff-mean_diff, 'Difference of mean SSH between data and model - offset', 'MADT [m]', plt.cm.bwr, plot_path, 'diff-offset_bwr.png')

# Diff - offset : change colormap
func_plot.polar_map(np.linspace(-0.5,0.5,50), Xa, Ya, diff-mean_diff, 'Difference of mean SSH between data and model - offset', 'MADT [m]', plt.cm.jet, plot_path, 'diff-offset_jet.png')

# Diff - offset / std_aviso
func_plot.polar_map(np.linspace(-5,14,50), Xa, Ya, (diff-mean_diff)/std_aviso, 'Difference of mean SSH between data and model - mean_diff / std AVISO', 'MADT [m]', plt.cm.jet, plot_path, 'diff_mean_diff_std_aviso.png')

# Diff / std_aviso
func_plot.polar_map(np.linspace(-5,14,50), Xa, Ya, diff/std_aviso, 'Difference of mean SSH between data and model / std AVISO', 'MADT [m]', plt.cm.jet, plot_path, 'diff_std_aviso.png')
