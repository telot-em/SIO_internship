import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap

def plot_map(X, Y, data, title, xlabel, ylabel, title_colorbar, path_save, title_save):
	'''
	Function for plotting a map.
	
	Args : 
	X, Y, data, title, xlabel, ylabel, title_colorbar, path_save, title_save
	'''

	plt.figure(figsize=(7,7))
	plt.contourf(X, Y, data)
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	c=plt.colorbar()
	c.set_label(title_colorbar)
	plt.savefig(path_save + title_save)
	plt.show()
	


def polar_map(cbar_levels, X, Y, data, title, title_colorbar, cmap, path_save, title_save):
	'''
	Function for plotting a polar projection map.
	
	Args : 
	cbar_levels, X, Y, data, title, title_colorbar, path_save, title_save
	'''
	
	# levels
	v       = cbar_levels

	# figure 
	fig     = plt.figure(figsize=(10,10))
	m       = Basemap(projection='ortho', lat_0=-90, lon_0=0)
	xm, ym  = m(X, Y)
	im      = m.contourf(xm, ym, data, levels=v, extend='both', cmap=plt.cm.jet)
	
	# colorbar
	cbar = m.colorbar(im, pad='20%')
	cbar.set_label(title_colorbar, fontsize=16)
	
	# structure
	m.fillcontinents(color='0.5', lake_color='0.5')
	m.drawparallels(np.arange(-90, -30, 20))#,labels=[True, False, True, True])
	m.drawmeridians(np.arange(0, 360, 30), labels=[1, 1, 1, 1]) 

	plt.title(title,fontsize=18, y=1.08)
	
	plt.show()
	if path_save is not None:
		plt.savefig(path_save + title_save, bbox_inches='tight')	
	
