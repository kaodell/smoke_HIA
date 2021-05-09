# regrid_population.py
# adapted from regrid.py code from Kelsey Bilsback
# January 2019
# Regrid satelite surface albedo data for RF calculations

# edited by Kate to use kriging grid and create US population

##################################
# user inputs
##################################
# grid to re-grid to
to_grid_fn = '/Users/kodell/Desktop/data_not4gdrive/kriging/kriging data/kriging_v2/krigedPM25_2006_v2_statfix_medbk.nc'
# population to regrid
pop_filename ='/Users/kodell/Desktop/data_not4gdrive/HIA_inputs/pop/gpw-v4-population-density-rev11_totpop_2pt5_min_nc/gpw_v4_population_density_rev11_2pt5_min.nc'
# out file and figure paths
fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/regrid_pop_check/'
out_file_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/population/'

##################################
# Import libraries
##################################
import numpy as np
import pylab as pl
import netCDF4
import datetime
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import interp
import matplotlib.colors as colors
#from area import area

def mplot():
    m = Basemap(llcrnrlon=-125,llcrnrlat=25,urcrnrlon=-65,urcrnrlat=50,
            projection='mill',lat_1=33,lat_2=45,lon_0=-95,resolution='l')
    # draw states
    m.drawstates()
    m.drawcoastlines()
    m.drawcountries()
    pl.box(on=None)
    return m

#############################################
# Read in, regrid, and combine netCDF file
#############################################
nc_fid = netCDF4.Dataset(to_grid_fn, 'r') # open file for reading
glon = nc_fid.variables['lon'][:].data
glat = nc_fid.variables['lat'][:].data
glon_we = nc_fid.variables['we_lon'][:].data
glat_we = nc_fid.variables['we_lat'][:].data
glon_ns = nc_fid.variables['ns_lon'][:].data
glat_ns = nc_fid.variables['ns_lat'][:].data
nc_fid.close()

nc_fid = netCDF4.Dataset(pop_filename, 'r') # open file to read
# data to regrid
pop_density_wm = nc_fid.variables['Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][:] # people per km2 
lon = nc_fid.variables['longitude'][:]
lat = nc_fid.variables['latitude'][:]
nc_fid.close()

pop_density_wm = pop_density_wm[2] # 2 = 2010
# replace fill value with zero
pop_density = np.ma.getdata(pop_density_wm)
mask = np.ma.getmask(pop_density_wm)
pop_density[mask] = 0.0

# cut this grid to the US
US_ind_latmax = np.max(np.where(lat > 15))
US_ind_latmin = np.min(np.where(lat < 55))
US_ind_lonmax = np.max(np.where(lon < -55))
US_ind_lonmin = np.min(np.where(lon > -140))

US_lon = np.copy(lon[US_ind_lonmin:US_ind_lonmax])
US_lat = np.copy(lat[US_ind_latmin:US_ind_latmax])
US_pop_density = np.copy(pop_density[ US_ind_latmin:US_ind_latmax,US_ind_lonmin:US_ind_lonmax])

print('Running grid interpolation')
# regrid pop data
# lat has to be increasing, so flip
US_lat_flip = np.flipud(US_lat)
US_pop_density_flip = np.flipud(US_pop_density)
pop_density_cu = interp(US_pop_density_flip, US_lon, US_lat_flip, glon, glat, order=1) 
# oreder=1 is a bi-linear interpolation

# hard-code grid area of 225 km
pop_cu = pop_density_cu * 225.0

#############################################
# write to file
#############################################
print('Writing to file')
# write file
nc_w_fid = netCDF4.Dataset(out_file_path+'rg_population_2010.nc', 'w', clobber=True,  format='NETCDF4')
nc_w_fid.description = 'Population from SEDAC in 2010'
nc_w_fid.history = 'Created' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# define file dimensions
nc_w_fid.createDimension('time', None) # unlimited dimension
nc_w_fid.createDimension('gridx', 189)
nc_w_fid.createDimension('gridy', 309)

lat_w = nc_w_fid.createVariable('lat', np.float32, ('gridx','gridy',))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('gridx','gridy',))
pop_w = nc_w_fid.createVariable('population', np.float32, ('gridx','gridy',))
pop_density_w = nc_w_fid.createVariable('population_density', np.float32, ('gridx','gridy',))
grid_area_w = nc_w_fid.createVariable('grid_area', np.float32, ('gridx','gridy'))

lon_w[:] = glon
lat_w[:] = glat
pop_w[:,:] = pop_cu
pop_density_w[:,:] = pop_density_cu
grid_area_w[:,:] = 225.0

nc_w_fid.close()

#############################################
# Check population plots
#############################################
pl.close('all')
'''
fig = pl.figure(figsize=(20, 20), dpi=120, facecolor='w', edgecolor='w')
pl.subplot(131)
# shift glon, glat to llcrn
m = mplot()
cs = m.pcolormesh(glon,glat, pop_cu,latlon=True,cmap='magma')
pl.title('Population counts')
cbar = pl.colorbar(orientation='horizontal')
cbar.set_label('Population per bin')
'''
# shift glon, glat to llcrn for pcolormesh in basemap
dx = 0.5*(glon[1:,1:]-glon[:-1,:-1])
dy = 0.5*(glat[1:,1:]-glat[:-1,:-1])

hfdx = np.vstack([dx[[0],:],dx])
fdx = np.hstack([hfdx[:,[0]],hfdx])

hfdy = np.vstack([dy[[0],:],dy])
fdy = np.hstack([hfdy[:,[0]],hfdy])

pglon = glon-fdx
pglat = glat-fdy

pl.figure()
pl.subplot(121)
m = mplot()
cs = m.pcolormesh(pglon,pglat, pop_density_cu,latlon=True,cmap='magma',vmax=50,norm=colors.LogNorm())
pl.title('Population density regrid')
cbar = pl.colorbar(orientation='horizontal')
cbar.set_label('Population per km2')

pl.subplot(122)
US_LON,US_LAT = np.meshgrid(US_lon,US_lat)
m = mplot()
cs = m.pcolormesh(US_LON,US_LAT, US_pop_density,latlon=True,cmap='magma',vmax=50,norm=colors.LogNorm())
pl.title('Population density original')
cbar = pl.colorbar(orientation='horizontal')
cbar.set_label('Population per km2')
pl.savefig(fig_path + 'population_regrid_check.png')
pl.show()

#%% different grid area tests and US pop calculations - not used in paper
'''
# grid with assigned country locations
id_grid_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/country_grid.npy'

# mask values outside the US
US_pop_cu = np.copy(pop_cu)
US_pop_density_cu = np.copy(pop_density_cu)

id_grid = np.load(id_grid_fn)
nonUS_inds = np.where(id_grid!='United States')
USinds = np.where(id_grid=='United States')

US_pop_cu[nonUS_inds] = np.nan
US_pop_density_cu[nonUS_inds] = np.nan

nonUSpop_cu = np.copy(pop_cu)
nonUSpop_cu[USinds] = np.nan

# this gives areas higher than 225 km? but if the WRF-Chem grid is 15x15km 
# not sure what is up with this; could have to do with the projection; corners being slightly off bc WRF is weird?

# calc area using code from Kelsey
# find grid corners
dx = 0.5*(glon[1:,1:]-glon[:-1,:-1])
dy = 0.5*(glat[1:,1:]-glat[:-1,:-1])

hfdx = np.vstack([dx[[0],:],dx])
fdx = np.hstack([hfdx[:,[0]],hfdx])

hfdy = np.vstack([dy[[0],:],dy])
fdy = np.hstack([hfdy[:,[0]],hfdy])

ll_crn_lat = glat-fdy
ll_crn_lon = glon-fdx

ur_crn_lat = glat+fdy
ur_crn_lon = glon+fdx

#grid_area_use = area_cu


area_cu = np.zeros([glat.shape[0],glat.shape[1]])
for i in range(glat.shape[0]):
        for j in range(glat.shape[1]):
                obj = {'type':'Polygon','coordinates':[[[glon[i,j]+fdx[i,j], glat[i,j]+fdy[i,j]],
		                                        [glon[i,j]+fdx[i,j], glat[i,j]-fdy[i,j]],
		                                        [glon[i,j]-fdx[i,j], glat[i,j]-fdy[i,j]],
		                                        [glon[i,j]-fdx[i,j], glat[i,j]+fdy[i,j]],
		                                        [glon[i,j]+fdx[i,j], glat[i,j]+fdy[i,j]]]]}
                area_cu[i,j] = np.float(area(obj)) * 1e-6 # area per km2
        print(i)

# calc area as dx x dy using midpoints of wrf grid
def haversine(lon0,lon1,lat0,lat1):    
    # from Dr. Will Lassman
    r = 6371000.#m                                                                                                                                                                                                                                                 
    lon0 = lon0*np.pi/180
    lon1 = lon1*np.pi/180
    lat0 = lat0*np.pi/180
    lat1 = lat1*np.pi/180
    return 2*r*np.arcsin(np.sqrt(np.sin((lat1 - lat0)/2.)**2 +\
		 np.cos(lat0)*np.cos(lat1)*np.sin((lon1 - lon0)/2.)**2))

        
gdx = haversine(glon_we[:,1:],glon_we[:,:-1],glat_we[:,1:],glat_we[:,:-1]).data
gdy = haversine(glon_ns[1:,:],glon_ns[:-1,:],glat_ns[1:,:],glat_ns[:-1,:]).data

gdx_km = gdx/1000.0
gdy_km = gdy/1000.0
grid_area = gdx_km*gdy_km

pop_cu_darea = pop_density_cu * grid_area
pop_cu_carea = pop_density_cu * area_cu

US_pop_cu_darea = np.copy(pop_cu_darea)
US_pop_cu_carea = np.copy(pop_cu_carea)

US_pop_cu_darea[nonUS_inds] = np.nan
US_pop_cu_carea[nonUS_inds] = np.nan


#############################################
# compare grid area methods and population differences
#############################################
fig = pl.figure()
pl.subplot(111)
m = mplot()
cs = m.pcolormesh(pglon,pglat, grid_area,latlon=True,cmap='magma')
pl.title('grid area')
cbar = pl.colorbar(orientation='horizontal')
cbar.set_label('area [km]')
pl.savefig('/Users/kodell/Desktop/grid_area_wrf.png')
'''
