import numpy as np
import netCDF4
import os
import pickle
import h5py

def pickle_save(name, path, data, verbose=True):
    if not os.path.exists(path):
        os.makedirs(path)
    full_name = (os.path.join(path,name+ '.npy'))


    with open(full_name, 'wb') as f2:
        pickle.dump(data, f2)
    if verbose:
        print('saved at : ',full_name)

load_path2='/data/SO12/runs/RUN_BLING_Dec2017/SO12_RUN/DIAGNOSTICS/'
load_path3='/data/soccom/GRID_12/'
folder = '/home/ebent/Octopus/Octopus-master/grid/'


f = h5py.File(load_path3 + 'grid.mat','r')
DYG = np.array(f.get('DYG'))
hFacW = np.array(f.get('hFacW'))
DRF = np.squeeze(np.array(f.get('DRF')))

#lon_min   = 900 
#lon_max   = 902 # 3240 si comme dans Traj_big_domain
#lat_min   = 512
#lat_max   = 514 # 1023 si comme dans Traj_big_domain
'''
nz,ny,nx = 104,1024,1801

fn1 = 'DYG.data'
fn2 = 'hFacW.data'
fn3 = 'DRF.data'

# 2D
DYG = np.fromfile(os.path.join(folder, fn1),'>f4')
DYG = np.reshape(DYG, [ny,nx])

# 3D
hFacW = np.fromfile(os.path.join(folder, fn2),'>f4')
hFacW = np.reshape(hFacW, [nz,ny,nx])

# 1D
DRF = np.fromfile(os.path.join(folder, fn3),'>f4')
'''


# This is hFacC for the SOUTHERN HEMISPHERE
file_h = h5py.File(load_path3 + 'grid.mat','r')

hFacC = file_h.get('hFacC')
hFacC = np.array(hFacC)
Xf = file_h.get('XC')
Xf = np.array(Xf)
Yf = file_h.get('YC')
Yf = np.array(Yf)

# Crop to the domain I want
lon_min   = 1320 # 110 deg ------- #1440 : same size as in Matlab
lon_max   = 3601 # 300 deg ------- #3241 : same size as in Matlab # 3240 : si comme dans Traj_big_domain
lat_min   = 0 
lat_max   = 1202 # -25 deg ------- #1024 : same size as in Matlab # 1023 si comme dans Traj_big_domain

hfacc = hFacC[:, lat_min:lat_max, lon_min:lon_max]
DYG = DYG[lat_min:lat_max, lon_min:lon_max]
hFacW = hFacW[:,lat_min:lat_max, lon_min:lon_max]
#DYG = DYG[lat_min:lat_max, lon_min:lon_max]
#hFacW = hFacW[:, lat_min:lat_max, lon_min:lon_max]

Uvel_list_names = ['so12_i0_year2006_5day_Uvel.nc', 
                   'so12_i0_year2007_5day_Uvel.nc', 
                   'so12_i0_year2008_5day_Uvel.nc', 
                   'so12_i0_year2009_5day_Uvel.nc', 
                   'so12_i0_year2010_5day_Uvel.nc', 
                   'so12_i0_year2011_5day_Uvel.nc']

# Calculate mean Uvel for k=0 (at surface)
mean_uvel = np.zeros((104,lat_max-lat_min,lon_max-lon_min))
tmp1 = np.zeros((lat_max-lat_min,lon_max-lon_min))
for name in Uvel_list_names:
	uvel = netCDF4.Dataset(load_path2 + str(name),'r')
	
	if name == 'so12_i0_year2006_5day_Uvel.nc':
		tmp1 = np.ma.mean(uvel.variables['Uvel'][19:,0,lat_min:lat_max, lon_min:lon_max], axis=0)
	else:
		tmp1 = tmp1 + np.ma.mean(uvel.variables['Uvel'][:,0,lat_min:lat_max, lon_min:lon_max], axis=0)
	mean_uvel[0] = tmp1/6

# Initialise tmp2 to calculate SF
tmp2 = mean_uvel[0]*DYG*hFacW[0,:,:]*DRF[0]

# Loop through all depths
for k in range(1,mean_uvel.shape[0]):
	print k
	for name in Uvel_list_names: # Calculate mean UVEL
    		print name
		uvel = netCDF4.Dataset(load_path2 + str(name),'r')
    		if name == 'so12_i0_year2006_5day_Uvel.nc':
        		tmp1 = np.ma.mean(uvel.variables['Uvel'][19:,k,lat_min:lat_max, lon_min:lon_max], axis=0)
			#print tmp1
    		else:
        		tmp1 = tmp1 + np.ma.mean(uvel.variables['Uvel'][:,k,lat_min:lat_max, lon_min:lon_max], axis=0)
			#print tmp1
	mean_uvel[k] = tmp1/6                                                                     
	tmp3 = mean_uvel[k]
	tmp3[np.where(tmp3==np.nan)] = 0
	mean_uvel[k] = tmp3

	tmp2 = tmp2 + tmp3*DYG*hFacW[k,:,:]*DRF[k] # integrate in z
        #print tmp2

pickle_save('tmp2_for_SF_3', '/data/ebent/', tmp2)

tmp2 = np.cumsum(tmp2, axis=0) # cumulative integral in y
tmp2[np.where(hfacc[0,...]==0)] = np.nan # masks 0s
SF = tmp2*1e-6 # converts from m3/s to Sv
pickle_save('SF_3', '/data/ebent/', SF)
pickle_save('mean_uvel_3', '/data/ebent', mean_uvel)


