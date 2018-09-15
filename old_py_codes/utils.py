def interp_etopo(tfile, lon_topo, lat_topo, z_topo, sub, lon_grid, lat_grid):
    '''
    Function to interpolate bathymetry grid to my data grid in order to then create a mask for topography >= 0 
    
    Args : tfile (file of topo/bathymetry), lon_topo, lat_topo, z_topo, sub (step used to select only part of 
    the data, lon_grid (lon of grid of our data on which we wish to interpolate topo), lat_grid (same).
    '''
    import numpy as np
    from netCDF4 import Dataset, num2date
    from scipy import interpolate

    etopo = Dataset(tfile, 'r')
    lats = etopo.variables[lat_topo][::sub] # selectes lats with a step of sub (to select all choose sub = 1)
    lons = etopo.variables[lon_topo][::sub]
    z = etopo.variables[z_topo][::sub,::sub]
    if ((lons >= -180.).all() and (lons <= 180.).all()): # .all() gets rid of Falses, only keeps Trues
        lon0 = np.where(lons>=0)[0][0] # to change lons from (-180, 180) to (0, 360)
        lons = lons + 180
        topo = z.copy()
        tmp = z[0,lon0:].shape[0]
        topo[:,:tmp] = z[:,lon0:]
        topo[:,tmp:] = z[:,:lon0]
    else: 
        topo = z.copy()

    grid = interpolate.interp2d(lons, lats, topo, kind='linear') # interpolates the topo to the grid of the data
    lon, lat = lon_grid, lat_grid
    etopo_interp = grid(lon, lat)
    return etopo_interp
