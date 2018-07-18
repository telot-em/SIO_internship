'''
Code written by Emma Bent on April 9, 2018

This code loads the data from AVISO (6 Apr 2006 to 5 Apr 2011).

Calculations made : 
	standard deviation of SSH (with land mask obtained thanks to ETOPO file)
	weighting of each cell by its area
	mean of SSH (with land mask as well)

Saves made : 
	std_aviso
	mean_aviso
	Xa (meshgrid of longitude)
	Ya (meshgrid of latitude)
	dsave_5 (SSH of AVISO selected every 5 days to create ice mask on 1/12 Model output that calculates SSH every 5 days (AVISO measures SSH every day))
	count_data (counts how many time steps of data are on each grid cell)
'''


# Import packages
import matplotlib.pyplot as plt
import numpy as np
import h5py
import utils
from netCDF4 import Dataset
from scipy import interpolate

load_path='/data/mmazloff/AVISO/'

# Change from .mat to numpy array
f = h5py.File(load_path + 'aviso_tot_MADT_sose_6Apr2006to5Apr2011.mat','r')

# Latitude and longitude variables
xsave = f.get('xsave')
xsave = np.array(xsave)
ysave = f.get('ysave')
ysave = np.array(ysave)

# Get rid of dimension 1 (go from 2D to 1D)
xsave = np.squeeze(xsave)
ysave = np.squeeze(ysave)

# Extract of data
tsave = f.get('tsave')
tsave = np.array(tsave)
tsave = np.squeeze(tsave)
dsave = f.get('dsave')
dsave = np.array(dsave)
dsave = np.moveaxis(dsave, 2, 0)
dsave = np.moveaxis(dsave, 2, 1)

print('Data has been loaded')

# Create grid for lat (Xa) and lon (Ya)
Xa, Ya = np.meshgrid(xsave, ysave)

#___________________________________________________________________________________#

# Calculate standard deviation of SSH
std_a = np.nanstd(dsave, axis=0)

#___________________________________________________________________________________#

# LAND MASK : interpolate topography on AVISO grid first

# ETOPO file
tfile = '/project_shared/ETOPO/ETOPO1_Ice_g_gmt4.grd'

# Interpolation of bathymetry on AVISO grid
etopo_interp = utils.interp_etopo(tfile=tfile, lon_topo='x', lat_topo='y', z_topo='z', sub=4, lon_grid=xsave, lat_grid=ysave)

# Applying a mask to the std deviation matrix
std_aviso = np.ma.masked_where(etopo_interp>=0, std_a)

#___________________________________________________________________________________#

# Taking into account the area of each cell which increases with latitude
cell_area = 0.25*0.25*np.cos(ysave*np.pi/180) # dx*dy*cos(lat), 0.25 = 1/4 of degree --> resolution of AVISO

dsave2 = np.moveaxis(dsave, 2, 1)
weight = dsave2*cell_area

weight_mean_aviso = np.nansum(weight, axis=0)/(cell_area*len(tsave))
weight_mean_aviso = weight_mean_aviso.T

# LAND MASK
weight_mean_aviso = np.ma.masked_where(etopo_interp>=0, weight_mean_aviso)

#___________________________________________________________________________________#

# Select part of dsave to make ICE MASK on model
dsave_5 = np.ma.masked_all((365, len(ysave), len(xsave)))
a = 0
b = 5
for c in range(365):
    dsave_5[c]= np.nanmean(dsave[a:b], axis=0)
    a+=5
    b+=5
dsave_5[-1,:,:]=dsave_5[-2,:,:]

#___________________________________________________________________________________#

# Count how many time steps of data are on each grid cell 
'''
count=np.zeros((len(ysave), len(xsave)))
for t in range(len(tsave)):
    for i in range(len(ysave)):
        for j in range(len(xsave)):
            if np.isnan(dsave[t,i,j]):
                count[i,j]=count[i,j]
            else:
                count[i,j]=count[i,j]+1
'''
#___________________________________________________________________________________#

# Create netCDF file to save
nc = Dataset('aviso_save.nc', 'w', format='NETCDF4')
nc.title = 'AVISO variables : std, mean, Xa, Ya, count_data'
from datetime import datetime
nc.date_created = datetime.now().isoformat()

nc.createDimension('lat', len(ysave))
nc.createDimension('lon', len(xsave))
nc.createDimension('time_model', 365)

std = nc.createVariable('std_aviso', 'f', ('lat','lon'))
std.long_name = 'Standard deviation of SSH of AVISO data'
std[:] = std_aviso

mean = nc.createVariable('mean_aviso', 'f', ('lat','lon'))
mean.long_name = 'Mean SSH of AVISO data'
mean[:] = weight_mean_aviso

X = nc.createVariable('Xa', 'f', ('lat','lon'))
X.long_name = 'Longitude meshgrid'
X[:] = Xa

Y = nc.createVariable('Ya', 'f', ('lat','lon'))
Y.long_name = 'Latitude meshgrid'
Y[:] = Ya


dsave_5 = nc.createVariable('dsave_5', 'f', ('time_model', 'lat', 'lon'))
dsave_5.long_name = 'AVISO selected every 5 time step to have time same size as model'
dsave_5[:] = dsave_5
'''
count_data = nc.createVariable('count', 'f', ('lat','lon'))
count_data.long_name = 'Count how many time steps of data on each grid cell'
count_data[:] = count
'''
nc.close()



'''

nc = Dataset('aviso_try.nc', 'w', format='NETCDF4')
nc.createDimension('lat', len(ysave))
nc.createDimension('lon', len(xsave))
std = nc.createVariable('std_aviso', '>f4', ('lat','lon'))
std[:] = std_aviso
X = nc.createVariable('Xa', 'f', ('lat','lon'))
Y = nc.createVariable('Ya', 'f', ('lat','lon'))
X[:] = Xa
Y[:] = Ya

nc.close()

f = Dataset('aviso_try.nc', 'r', format='NETCDF4')
std_avisooo = f.variables['std_aviso'][:]
'''
 
'''
# Extract a 'chunk' of whole data
data_chunk = f['dsave'][:,:,:10]

chunk_mean = np.nanmean(data_chunk, axis=2)
std_aviso = np.nanstd(data_chunk, axis=2)
'''

#------------------------PLOTS---------------------------
'''
## Plots for period 2006_2011

# Plot the mean SSH
my_functions.plot_map(X, Y, mean_ssh, 'Mean SSH for 2006_2011', 'Longitude', 'Latitude', 'MADT [m]', plot_path_2006_2011, 'mean_ssh.png')

# Plot the mean SSH for chunk
my_functions.plot_map(X, Y, chunk_mean, 'Mean SSH for chunk', 'Longitude', 'Latitude', 'MADT [m]', plot_path_2006_2011, 'ssh_chunk.png')
'''
#-------------------POLAR PROJ---------------------------

# Plot the mean SSH - POLAR PROJ
#my_functions.polar_map(np.linspace(-2,2,50), Xa, Ya, mean_ssh_aviso, 'Mean SSH for 2006_2011(AVISO data)', 'MADT [m]', plot_path_2006_2011, 'aviso_polar_mean_ssh.png') 

# Plot the mean SSH of chunk - POLAR PROJ
#my_functions.polar_map(np.linspace(-2,2,50), Xa, Ya, chunk_mean, 'Mean SSH for chunk (AVISO data)', 'MADT [m]', plot_path_2006_2011, 'polar_chunk.png')

# Plot standard deviation
#my_functions.polar_map(np.linspace(0, 0.52, 50), Xa, Ya, std_aviso, 'Standard deviation for 2006_2011 (AVISO data)', 'Standard deviation [m]', plot_path_2006_2011, 'aviso_std.png')


