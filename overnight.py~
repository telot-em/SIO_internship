import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, pi
import h5py
from scipy.io import loadmat
from mpl_toolkits.basemap import Basemap
import netCDF4
from scipy import interpolate
import os
import pickle

def pickle_save(name, path, data, verbose=True):
    if not os.path.exists(path):
        os.makedirs(path)
    full_name = (os.path.join(path,name+ '.npy'))


    with open(full_name, 'wb') as f2:
        pickle.dump(data, f2)
    if verbose:
        print('saved at : ',full_name)

def pickle_load(name, path, verbose=True):
    #if not os.path.exists(path):
    #    os.makedirs(path)
    full_name= (os.path.join(path,name+ '.npy'))

    with open(full_name, 'r') as f:
        data=pickle.load(f)

    if verbose:
        print('loaded from : ',full_name)
    return data

plot_path_1993_2017='/home/ebent/plots/1993_2017/'
plot_path_2006_2011='/home/ebent/plots/2006_2011/'
plot_path_jup='/home/ebent/plots/2006_2011/jup2/'
load_path='/data/mmazloff/AVISO/'
load_path2='/data/SO12/runs/RUN_BLING_Dec2017/SO12_RUN/DIAGNOSTICS/'
load_path3='/data/soccom/GRID_12/'

file1 = netCDF4.Dataset(load_path2+'so12_i0_year2006_5day_Salt.nc','r')
file2 = netCDF4.Dataset(load_path2+'so12_i0_year2007_5day_Salt.nc','r')
file3 = netCDF4.Dataset(load_path2+'so12_i0_year2008_5day_Salt.nc','r')
file4 = netCDF4.Dataset(load_path2+'so12_i0_year2009_5day_Salt.nc','r')
file5 = netCDF4.Dataset(load_path2+'so12_i0_year2010_5day_Salt.nc','r')
file6 = netCDF4.Dataset(load_path2+'so12_i0_year2011_5day_Salt.nc','r')
'''
file1 = netCDF4.Dataset(load_path2+'so12_i0_year2006_5day_Theta.nc','r')
file2 = netCDF4.Dataset(load_path2+'so12_i0_year2007_5day_Theta.nc','r')
file3 = netCDF4.Dataset(load_path2+'so12_i0_year2008_5day_Theta.nc','r')
file4 = netCDF4.Dataset(load_path2+'so12_i0_year2009_5day_Theta.nc','r')
file5 = netCDF4.Dataset(load_path2+'so12_i0_year2010_5day_Theta.nc','r')
file6 = netCDF4.Dataset(load_path2+'so12_i0_year2011_5day_Theta.nc','r')
'''
#choice of mercator box, so all region

lon_180W = 2160
lon_150W = 2520
lat_min = 0
lat_max = 1170
'''
lon_min = 1950
lon_max = 2520
lat_min = 0
lat_max = 541
'''
'''
lat = file1.variables['lat'][lat_min:lat_max]
lon = file1.variables['lon'][:]
time = file1.variables['time'][:]
depth = file1.variables['depth'][:]

mean_Theta_150W= (np.ma.mean(file1.variables['Theta'][19:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file2.variables['Theta'][:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file3.variables['Theta'][:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file4.variables['Theta'][:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file5.variables['Theta'][:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file6.variables['Theta'][:,:,lat_min:lat_max, lon_150W], axis=0))/6

pickle_save('mean_Theta_150W', '/data/ebent', mean_Theta_150W)

mean_Theta_180W= (np.ma.mean(file1.variables['Theta'][19:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file2.variables['Theta'][:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file3.variables['Theta'][:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file4.variables['Theta'][:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file5.variables['Theta'][:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file6.variables['Theta'][:,:,lat_min:lat_max, lon_180W], axis=0))/6

pickle_save('mean_Theta_180W', '/data/ebent', mean_Theta_180W)

'''
mean_Salt_150W= (np.ma.mean(file1.variables['Salt'][19:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file2.variables['Salt'][:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file3.variables['Salt'][:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file4.variables['Salt'][:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file5.variables['Salt'][:,:,lat_min:lat_max, lon_150W], axis=0) + \
np.ma.mean(file6.variables['Salt'][:,:,lat_min:lat_max, lon_150W], axis=0))/6

pickle_save('mean_Salt_150W', '/data/ebent', mean_Salt_150W)

mean_Salt_180W= (np.ma.mean(file1.variables['Salt'][19:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file2.variables['Salt'][:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file3.variables['Salt'][:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file4.variables['Salt'][:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file5.variables['Salt'][:,:,lat_min:lat_max, lon_180W], axis=0) + \
np.ma.mean(file6.variables['Salt'][:,:,lat_min:lat_max, lon_180W], axis=0))/6

pickle_save('mean_Salt_180W', '/data/ebent', mean_Salt_180W)

