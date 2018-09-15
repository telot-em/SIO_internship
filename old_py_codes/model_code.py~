'''
Code written by Emma Bent on April 9, 2018

This code loads the 1/12 Model output for SSH (6 Apr 2006 to 5 Apr 2011).

Calculations made : 
	interpolation of the model output to the AVISO grid for same period
	ice mask (obtained with AVISO) on the model output
        standard deviation of SSH (with land mask obtained thanks to hFacC file)
        weighting of each cell by its area + ice mask 
        mean of SSH (with land mask as well)
	

Saves made : 
        std_model 
        mean_model
'''


# Import packages
import matplotlib.pyplot as plt 
import numpy as np
import h5py
import utils 
from netCDF4 import Dataset
from scipy import interpolate

load_path='/data/SO12/runs/RUN_BLING_Dec2017/SO12_RUN/DIAGNOSTICS/'
load_path2='/data/mmazloff/AVISO/'
load_path3='/data/soccom/GRID_12/'

# Read data
f = Dataset(load_path+'SO12_Dec17_2005toAp2011_5day_SSH.nc','r')

# Latitude and longitude variables
XC=f.variables['XC'][:]
YC=f.variables['YC'][:]

# Extract of data
time=f.variables['time'][91:-1]
ssh=f.variables['ETAN'][91:-1]
area=f.variables['rA'][:]

print('Data has been loaded')

# Create grid for lat (Xm) and lon (Ym)
Xm, Ym = np.meshgrid(XC, YC) 

# Load variables saved from the AVISO data
g = Dataset('aviso_test2.nc', 'r', format='NETCDF4')

dsave_5=g.variables['dsave_5'][:]
xsave=g.variables['Xa'][0,:]
ysave=g.variables['Xa'][:,0]

#___________________________________________________________________________________#

# Interpolation of the whole model output on AVISO grid
ssh_interp = np.ma.masked_all((len(time), len(ysave), len(xsave)))
for t in range(len(time)):
    gridd = interpolate.interp2d(XC, YC, ssh[t], kind='linear') # XC and YC are lat and lon vectors of model
    ssh_interp[t] = gridd(xsave, ysave)

# ICE MASK

ssh_ice_mask = np.ma.masked_all((len(time), len(ysave), len(xsave)))
for l in range(len(time)):
    ssh_ice_mask[l] = np.ma.masked_where(np.isnan(dsave_5[l]), ssh_interp[l])

#___________________________________________________________________________________#

# Calculate standard deviation
std_ice = np.ma.std(ssh_ice_mask, axis=0)

#___________________________________________________________________________________#

# LAND MASK

file = h5py.File(load_path3 + 'grid.mat','r')

hFacC = file.get('hFacC')
hFacC = np.array(hFacC)
Xf = file.get('XC')
Xf = np.array(Xf)
Yf = file.get('YC')
Yf = np.array(Yf)
RC = file.get('RC')
RC = np.array(RC)
z = np.squeeze(RC)

x = Xf[0,:]
y = Yf[:,0]

# hFacC[0,:,:] corresponds to the layer at the surface (z=-1), we interpolate it to the AVISO grid
griD = interpolate.interp2d(x, y, hFacC[0,:,:], kind='linear') # XC and YC are lat and lon vectors of model
hFacC_interp = griD(xsave, ysave)

# We apply a land mask to the std of the model 
std_model = np.ma.masked_where(hFacC_interp==0, std_ice)

#___________________________________________________________________________________#

# Taking into account the area of each cell which increases with latitude
weight_ssh = ssh*area

# Interpolate area to AVISO grid
gridd = interpolate.interp2d(XC, YC, area, kind='linear')
area_interp = gridd(xsave, ysave)

# Interpolate weight_ssh to AVISO grid
weight_ssh_interp = np.ma.masked_all((len(time), len(ysave), len(xsave)))
for t in range(len(time)):
    gridd = interpolate.interp2d(XC, YC, weight_ssh[t], kind='linear') # XC and YC are lat and lon vectors of model
    weight_ssh_interp[t] = gridd(xsave, ysave)

# ICE MASK
weight_ice_mask = np.ma.masked_all((len(time), len(ysave), len(xsave)))
for l in range(len(time)):
    weight_ice_mask[l] = np.ma.masked_where(np.isnan(dsave_5[l]), weight_ssh_interp[l])

# Calculate the mean
weight_mean_model = np.ma.sum(weight_ssh_interp, axis=0)/(area_interp*len(time))

# LAND MASK
weight_mean_model = np.ma.masked_where(hFacC_interp==0, weight_mean_model)

#___________________________________________________________________________________#

# Create netCDF file to save
nc = Dataset('model_save.nc', 'w', format='NETCDF4')
nc.title = 'Model  variables : std, mean'
from datetime import datetime
nc.date_created = datetime.now().isoformat()

nc.createDimension('lat', len(ysave))
nc.createDimension('lon', len(xsave))

std = nc.createVariable('std_model', 'f', ('lat','lon'))
std.long_name = 'Standard deviation of SSH of Model'
std[:] = std_model

mean = nc.createVariable('mean_model', 'f', ('lat','lon'))
mean.long_name = 'Mean SSH of Model'
mean[:] = weight_mean_model

nc.close()



'''
#------------------------PLOTS---------------------------

# Plot the mean SSH - POLAR PROJ
#my_functions.polar_map(np.linspace(-2,2,50), Xm, Ym, mean_ssh_model, 'Mean SSH for 2006_2011 (model)', 'MADT [m]', plot_path_2006_2011, 'model_polar_mean_ssh.png') 

# Plot the mean SSH of chunk - POLAR PROJ
#my_functions.polar_map(np.linspace(-2,2,50), Xm, Ym, chunk_mean, 'Mean SSH for chunk (model)', 'MADT [m]', plot_path_2006_2011, 'model_polar_chunk.png')

# Plot standard deviation
my_functions.polar_map(np.linspace(0, 0.52, 50), Xm, Ym, std_model, 'Standard deviation for 2006_2011 (model)', 'Standard deviation [m]', plot_path_2006_2011, 'model_std.png')


'''

