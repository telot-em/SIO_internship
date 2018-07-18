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

plot_path_1993_2017='/home/ebent/plots/1993_2017/'
plot_path_2006_2011='/home/ebent/plots/2006_2011/'
plot_path_jup='/home/ebent/plots/2006_2011/jup2/'
load_path='/data/mmazloff/AVISO/'
load_path2='/data/SO12/runs/RUN_BLING_Dec2017/SO12_RUN/DIAGNOSTICS/'
load_path3='/data/soccom/GRID_12/'
file1 = netCDF4.Dataset(load_path2+'so12_i0_year2006_5day_Theta.nc','r')
file2 = netCDF4.Dataset(load_path2+'so12_i0_year2007_5day_Theta.nc','r')
file3 = netCDF4.Dataset(load_path2+'so12_i0_year2008_5day_Theta.nc','r')
file4 = netCDF4.Dataset(load_path2+'so12_i0_year2009_5day_Theta.nc','r')
file5 = netCDF4.Dataset(load_path2+'so12_i0_year2010_5day_Theta.nc','r')
file6 = netCDF4.Dataset(load_path2+'so12_i0_year2011_5day_Theta.nc','r')
# Select a specific region
lon_min = 1950
lon_max = 2520
lat_min = 0
lat_max = 541
lat = file1.variables['lat'][lat_min:lat_max]
lon = file1.variables['lon'][lon_min:lon_max]
Lon, Lat = np.meshgrid(lon,lat)
depth = file1.variables['depth'][:]
surf = depth[0]
depth_100 = depth[24]
depth_200 = depth[39]
print(surf, depth_100, depth_200)
