'''
Code written by Emma Bent on April 9, 2018

This code loads mean SSH from AVISO data and 1/12 Model output for the period 6 Apr 2006 to 5 Apr 2011 and calculates their difference and the mean of that difference to know the offset between the data and the model

Saves made :
	diff
	mean_diff
'''

from netCDF4 import Dataset
import numpy as np

f = Dataset('aviso_save.nc', 'r', format='NETCDF4')
weight_mean_aviso = f.variables['mean_aviso'][:]
ysave = f.variables['Ya']
xsave = f.variables['Xa']
ysave = ysave[:,0]
xsave = xsave[0,:]

f2 = Dataset('model_save.nc', 'r', format='NETCDF4')
weight_mean_model = f2.variables['mean_model'][:]


# Difference between mean SSH of data and model
diff_weight = weight_mean_aviso-weight_mean_model

# Calculate the offset : mean_diff
mean_diff_weight = np.nanmean(diff_weight)

# Create netCDF file to save
nc = Dataset('diff_save.nc', 'w', format='NETCDF4')

nc.title = 'Difference between mean SSH of AVISO data and 1/12 Model'

from datetime import datetime
nc.date_created = datetime.now().isoformat()

nc.createDimension('lat', len(ysave))
nc.createDimension('lon', len(xsave))

diff = nc.createVariable('diff_weight', 'f', ('lat','lon'))
diff.long_name = 'Difference between mean SSH of AVISO data and 1/12 Model' 
diff = diff_weight

mean_diff = nc.createVariable('mean_diff_weight', 'f', ('lat','lon'))
mean_diff.long_name = 'Mean of difference between mean SSH of data and model = offset'
mean_diff = mean_diff_weight

nc.close()
