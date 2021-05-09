# Calculate mortalities
# Bilsback 2020
# edited by Kate to ouput maps of chronic smoke PM and population
import numpy as np
import netCDF4
import datetime
import matplotlib as mplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors
import cartopy.io.shapereader as shpreader
from cartopy.feature import NaturalEarthFeature
from wrf import (getvar, get_cartopy, cartopy_xlim,cartopy_ylim)
#%% user inputs
# file locations and names
kPM_fp = '/Users/kodell/Desktop/data_not4gdrive/HIA_inputs/PM/'
kPM_desc = '_v2_medbk_final' # which kPM file to load
pop_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/population/rg_population_2010.nc'
state_grid_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/state_grid.npy'

# out file paths
out_fn_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/mortalities/'
out_fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/gridded_mort_maps/'
desc = kPM_desc+'_2010pop'

# years to load
years = np.arange(2006,2019)

#%% user-defined functions
def plot_background(ax):
    ax.set_extent([235., 290., 20., 55.])
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    return ax

def mk_map(ax):
    ax.patch.set_visible(False)
    shapename = 'admin_1_states_provinces_lakes_shp'
    states_shp = shpreader.natural_earth(resolution='110m',category='cultural', name=shapename)
    facecolor = (0, 0, 0, 0)
    edgecolor = 'black'
    for state in shpreader.Reader(states_shp).geometries():
        ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)

def calc_mort(disease, theta, se_theta, alpha, mu, pi, base,
              totalPM, nosmokePM,sf_avg_pm):

    print('Calc mortalities for ', disease)
    thetas = [theta - 2 * se_theta, theta, theta + 2 * se_theta]

    z_total = np.where(totalPM > 2.4, totalPM - 2.4, 0)
    z_nosmoke = np.where(nosmokePM > 2.4, nosmokePM - 2.4, 0)
    
    Gamma_total = np.log(1 + (z_total / alpha)) / (1 + np.exp((mu - z_total) / (pi)))
    Gamma_nosmoke = np.log(1 + (z_nosmoke / alpha)) / (1 + np.exp((mu - z_nosmoke) / (pi)))

    # Calculate hazard ratio
    HR_total = np.exp(np.array(thetas)[:,None,None] * Gamma_total)
    HR_nosmoke = np.exp(np.array(thetas)[:,None,None] * Gamma_nosmoke)

    #Mortalities
    M_total = base * population[:,:] * (1 - (1 / HR_total)) * 1e-5 # last number adjusts for baseline mortality units
    M_nosmoke = base * population[:,:] * (1 - (1 / HR_nosmoke)) * 1e-5

    dM_smoke = M_total - M_nosmoke

    #attributable fraction
    dM_af = sf_avg_pm*M_total

    return(dM_smoke, dM_af, M_total, M_nosmoke)
#%% load and prep PM data
print('Read in pm data')
#Get pm data and calc average
kPM_fid = netCDF4.Dataset(kPM_fp + 'krigedPM25_06-18'+kPM_desc+'.nc')
allday_smokePM = kPM_fid['Smoke PM25'][:].data
allday_bkPM = kPM_fid['Background PM25'][:].data
allday_totPM = kPM_fid['PM25'][:].data
HMS_smoke = kPM_fid['HMS Smoke'][:].data
lon = kPM_fid['lon'][:].data
lat = kPM_fid['lat'][:].data
kPM_fid.close()

# make non smoke array where on non-HMS days = total PM and on smoke days = bkPM
allday_nosmokePM = np.where(HMS_smoke==0,allday_totPM,allday_bkPM)

# the sum of smoke and nosmoke should equal the total PM, check
summ = allday_nosmokePM + allday_smokePM

diff = allday_totPM - summ
max_diff = np.nanmax(diff)
min_diff = np.nanmin(diff)
print('these numbers should be essentially zero: ',max_diff,min_diff)

# calc PM averages 
pm_total = np.nanmean(allday_totPM,axis=(0,1))
pm_nosmoke = np.nanmean(allday_nosmokePM,axis=(0,1))
pm_smoke = np.nanmean(allday_smokePM,axis=(0,1))
# calc smoke fraction average. Do this way vs avg(daily smoke/daily total) to 
# avoid issues with small #/small #
sf_avg = pm_smoke/pm_total

# load population
nc_pop = netCDF4.Dataset(pop_file, 'r')
population = nc_pop.variables['population'][:]
area = nc_pop.variables['grid_area'][:]
nc_pop.close()

# load state grid to remove values outside the US
state_grid = np.load(state_grid_fn)
non_US_inds = np.where(state_grid=='NA')
pm_total[non_US_inds] = np.nan
pm_nosmoke[non_US_inds] = np.nan
sf_avg[non_US_inds] = np.nan
population[non_US_inds] = np.nan
#%% calculate mortalities
print('calculating mortalities')
# inputs in function are: (disease, theta, se_theta, alpha, mu, pi, base, totalPM, nosmokePM)
# baseline values are from here: http://ghdx.healthdata.org/gbd-results-tool (2019 analysis)
# on webpage select: location: Location: US, Year: 2010 Metric: Rate, Age: All ages, Measure: Death, Sex: Both

# other function values are from the Burnett et al., 2018 supplementary file: GEMM calculator

# 5 causes
# COPD - Chronic Obstructive Pulmany Disease
[dM_smoke_copd, dM_af_copd, M_total_copd, M_nosmoke_copd] = calc_mort('COPD', 0.2510, 0.06762, 6.5, 2.5, 32.0, 51.27, pm_total, pm_nosmoke,sf_avg)
# IHD - Ischemic Heart Disease
[dM_smoke_ihd, dM_af_ihd, M_total_ihd, M_nosmoke_ihd] = calc_mort('IHD', 0.2969, 0.01787, 1.9, 12.0, 40.2, 158.63, pm_total, pm_nosmoke,sf_avg)
# Lung Cancer
[dM_smoke_lc, dM_af_lc, M_total_lc, M_nosmoke_lc] = calc_mort('LC', 0.2942, 0.06147, 6.2, 9.3, 29.8, 58.73, pm_total, pm_nosmoke,sf_avg)
# Stroke
[dM_smoke_stroke, dM_af_stroke, M_total_stroke, M_nosmoke_stroke] = calc_mort('STROKE', 0.2720, 0.07697, 6.2, 16.7, 23.7, 51.46, pm_total, pm_nosmoke,sf_avg)
# LRI - Lower Respiratory Infection
[dM_smoke_lri, dM_af_lri, M_total_lri, M_nosmoke_lri] = calc_mort('LRI',  0.4468, 0.11735, 6.4, 5.7, 8.4, 22.21, pm_total, pm_nosmoke,sf_avg)

# all-cause (non-communicalbe diseases + LRI)
[dM_smoke_allcause, dM_af_allcause, M_total_allcause, M_nosmoke_allcause] = calc_mort('ALL_CAUSE', 0.1430, 0.01807, 1.6, 15.5, 36.8, 732.93, pm_total, pm_nosmoke,sf_avg)

#%% calculate totals
print('Calculate total deaths')

dM_total_smoke = dM_smoke_copd + dM_smoke_ihd + dM_smoke_lc + dM_smoke_stroke + dM_smoke_lri
print('deaths prevented if no smoke, 5 leading causes')
print(np.nansum(dM_total_smoke,axis=(1,2)))

dM_total_af = dM_af_copd + dM_af_ihd + dM_af_lc + dM_af_stroke + dM_af_lri
print('deaths currently attributable to smoke, 5 leading causes')
print(np.nansum(dM_total_af,axis=(1,2)))

M_total_total = M_total_copd + M_total_ihd + M_total_lc + M_total_stroke + M_total_lri
M_nosmoke_total = M_nosmoke_copd + M_nosmoke_ihd + M_nosmoke_lc + M_nosmoke_stroke + M_nosmoke_lri

# compare to all-cause
print('deaths prevented if no smoke, all-cause')
print(np.nansum(dM_smoke_allcause,axis=(1,2)))

dM_total_af = dM_af_copd + dM_af_ihd + dM_af_lc + dM_af_stroke + dM_af_lri
print('deaths currently attributable to smoke, all-cause')
print(np.nansum(dM_af_allcause,axis=(1,2)))


#%% save to file
print('saving to file')
nc_w_fid = netCDF4.Dataset(out_fn_path+'wilfire_smoke_premature_mort_gemm_'+desc+'.nc', 'w', clobber=True,  format='NETCDF4')
nc_w_fid.description = 'wildfire smoke premature mortalities'
nc_w_fid.history = 'Created ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# define file dimensions
nc_w_fid.createDimension('lat', 189)
nc_w_fid.createDimension('lon', 309)
nc_w_fid.createDimension('stats', 3)

lat_w = nc_w_fid.createVariable('lat', np.float32, ('lat','lon'))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('lat','lon'))

# 5 leading causes
smoke_mort_w = nc_w_fid.createVariable('smoke_mort_averted', np.float32, ('stats', 'lat','lon',))
af_mort_w = nc_w_fid.createVariable('smoke_mort_att_fraction', np.float32, ('stats', 'lat','lon',))

smoke_mort_density_w = nc_w_fid.createVariable('smoke_mort_averted_density', np.float32, ('stats', 'lat','lon',))
af_mort_density_w = nc_w_fid.createVariable('smoke_mort_att_fraction_density', np.float32, ('stats', 'lat','lon',))

total_mort_w = nc_w_fid.createVariable('total_mort', np.float32, ('stats', 'lat','lon',))
nosmoke_mort_w = nc_w_fid.createVariable('nosmoke_mort', np.float32, ('stats', 'lat','lon',))

# all cause (non-communicable + lri)
smoke_mort_ac_w = nc_w_fid.createVariable('smoke_mort_averted_ac', np.float32, ('stats', 'lat','lon',))
af_mort_ac_w = nc_w_fid.createVariable('smoke_mort_att_fraction_ac', np.float32, ('stats', 'lat','lon',))

smoke_mort_density_ac_w = nc_w_fid.createVariable('smoke_mort_averted_density_ac', np.float32, ('stats', 'lat','lon',))
af_mort_density_ac_w = nc_w_fid.createVariable('smoke_mort_att_fraction_density_ac', np.float32, ('stats', 'lat','lon',))

total_mort_ac_w = nc_w_fid.createVariable('total_mort_ac', np.float32, ('stats', 'lat','lon',))
nosmoke_mort_ac_w = nc_w_fid.createVariable('nosmoke_mort_ac', np.float32, ('stats', 'lat','lon',))


lon_w[:] = lon
lat_w[:] = lat

# 5 leading causes
smoke_mort_w[:,:,:] = dM_total_smoke
af_mort_w[:,:,:] = dM_total_af

smoke_mort_density_w[:,:,:] = dM_total_smoke / area
af_mort_density_w[:,:,:] = dM_total_af / area

total_mort_w[:,:,:] = M_total_total
nosmoke_mort_w[:,:,:] = M_nosmoke_total

# all cause
smoke_mort_ac_w[:,:,:] = dM_smoke_allcause
af_mort_ac_w[:,:,:] = dM_af_allcause

smoke_mort_density_ac_w[:,:,:] = dM_smoke_allcause / area
af_mort_density_ac_w[:,:,:] = dM_af_allcause / area

total_mort_ac_w[:,:,:] = M_total_allcause
nosmoke_mort_ac_w[:,:,:] = M_nosmoke_allcause

nc_w_fid.close()

#%% calculate DALYs and save
print('Calculate DALYs')
# calc all-cause (non-communicalbe diseases + LRI)
[dD_smoke_allcause, dD_af_allcause, D_total_allcause, D_nosmoke_allcause] = calc_mort('ALL_CAUSE', 0.1430, 0.01807, 1.6, 15.5, 36.8, 26688.87, pm_total, pm_nosmoke,sf_avg)
dD_total_smoke = dD_smoke_allcause
dD_total_af = dD_af_allcause
D_total_total = D_total_allcause
D_nosmoke_total = D_nosmoke_allcause

nc_w_fid = netCDF4.Dataset(out_fn_path+'wilfire_smoke_PM_DALYs_gemm'+desc+'.nc', 'w', clobber=True,  format='NETCDF4')
nc_w_fid.description = 'wildfire smoke DALYs'
nc_w_fid.history = 'Created ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# define file dimensions
nc_w_fid.createDimension('lat', 189)
nc_w_fid.createDimension('lon', 309)
nc_w_fid.createDimension('stats', 3)

lat_w = nc_w_fid.createVariable('lat', np.float32, ('lat','lon'))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('lat','lon'))

smoke_DALY_w = nc_w_fid.createVariable('smoke_DALY_averted', np.float32, ('stats', 'lat','lon',))
af_DALY_w = nc_w_fid.createVariable('smoke_DALY_fraction', np.float32, ('stats', 'lat','lon',))

smoke_DALY_density_w = nc_w_fid.createVariable('smoke_DALY_averted_density', np.float32, ('stats', 'lat','lon',))
af_DALY_density_w = nc_w_fid.createVariable('smoke_DALY_att_fraction_density', np.float32, ('stats', 'lat','lon',))

total_DALY_w = nc_w_fid.createVariable('total_DALY', np.float32, ('stats', 'lat','lon',))
nosmoke_DALY_w = nc_w_fid.createVariable('nosmoke_DALY', np.float32, ('stats', 'lat','lon',))

lon_w[:] = lon
lat_w[:] = lat

smoke_DALY_w[:,:,:] = dD_total_smoke
af_DALY_w[:,:,:] = dD_total_af

smoke_DALY_density_w[:,:,:] = dD_total_smoke / area
af_DALY_density_w[:,:,:] = dD_total_af / area

total_DALY_w[:,:,:] = D_total_total
nosmoke_DALY_w[:,:,:] = D_nosmoke_total

nc_w_fid.close()

#%% make plots
# mortalities
# mask values outside the US, so we can get sum of US population
plot_Mtotal = (M_total_allcause[1,:,:])/area
plot_nsmokeMtotal = (M_nosmoke_allcause[1,:,:])/area
plot_Mtotal_smoke = (dM_smoke_allcause[1,:,:])/area
plot_afMsmoke = (dM_af_allcause[1,:,:])/area

# calc total US deaths for each so we can put on the map
US_Mtotal = np.nansum(M_total_allcause,axis=(1,2))
US_Mnosmoke = np.nansum(M_nosmoke_allcause,axis=(1,2))
US_Msmoke = np.nansum(dM_smoke_allcause,axis=(1,2))
US_Msmoke_af = np.nansum(dM_af_allcause,axis=(1,2))
# for 5 leading causes
US_Mtotal_5l = np.nansum(M_total_total,axis=(1,2))
US_Mnosmoke_5l = np.nansum(M_nosmoke_total,axis=(1,2))
US_Msmoke_5l = np.nansum(dM_total_smoke,axis=(1,2))
US_Msmoke_af_5l = np.nansum(dM_total_af,axis=(1,2))

total_str = str(round(US_Mtotal[0])) + '-' + str(round(US_Mtotal[-1]))
nosmoke_str = str(round(US_Mnosmoke[0])) + '-' + str(round(US_Mnosmoke[-1]))
smoke_str = str(round(US_Msmoke[0])) + '-' + str(round(US_Msmoke[-1]))
af_str = str(round(US_Msmoke_af[0])) + '-' + str(round(US_Msmoke_af[-1]))

ncfile = netCDF4.Dataset('/Users/kodell/Desktop/data_not4gdrive/wrfout_d01_2016-07-22_00_00_00')
slp = getvar(ncfile,"slp")
cart_proj = get_cartopy(slp)

fig, axarr = plt.subplots(nrows=2, ncols=2, subplot_kw={'projection': cart_proj})

axlist = axarr.flatten()
for ax in axlist:
    mk_map(ax)
cs=axlist[0].pcolor(lon,lat,plot_Mtotal, transform=ccrs.PlateCarree(),
        norm=colors.LogNorm(vmin=0.0001, vmax=10),shading='nearest')
axlist[0].set_title('total PM mortalities '+total_str, fontsize=16)
axlist[0].outline_patch.set_edgecolor('white')
cax,kw = mplt.colorbar.make_axes(axlist[0],location='bottom',pad=0.05,shrink=0.65)
cbar=fig.colorbar(cs,cax=cax,**kw)

cs=axlist[1].pcolor(lon,lat, plot_nsmokeMtotal, transform=ccrs.PlateCarree(),
                norm=colors.LogNorm(vmin=0.0001, vmax=10),shading='nearest')
axlist[1].set_title('total no smoke mortalties '+nosmoke_str, fontsize=16)
axlist[1].outline_patch.set_edgecolor('white')
cax,kw = mplt.colorbar.make_axes(axlist[1],location='bottom',pad=0.05,shrink=0.65)
cbar=fig.colorbar(cs,cax=cax,**kw)

cs=axlist[2].pcolor(lon,lat, plot_Mtotal_smoke, transform=ccrs.PlateCarree(),
                norm=colors.LogNorm(vmin=0.0001, vmax=1),shading='nearest')
axlist[2].set_title('deaths averted if no smoke '+smoke_str, fontsize=16)
axlist[2].outline_patch.set_edgecolor('white')
cax,kw = mplt.colorbar.make_axes(axlist[2],location='bottom',pad=0.05,shrink=0.65)
cbar=fig.colorbar(cs,cax=cax,**kw)

cs=axlist[3].pcolor(lon,lat, plot_afMsmoke, transform=ccrs.PlateCarree(),
                norm=colors.LogNorm(vmin=0.0001, vmax=1),shading='nearest')
axlist[3].set_title('deaths currently attributable to smoke '+af_str, fontsize=16)
axlist[3].outline_patch.set_edgecolor('white')
cax,kw = mplt.colorbar.make_axes(axlist[3],location='bottom',pad=0.05,shrink=0.65)
cbar=fig.colorbar(cs,cax=cax,**kw)

plt.savefig(out_fig_path + 'US_smoke_mort_'+desc+'.png')
plt.show()

# DALYs
plot_Dtotal = (D_total_total[1,:,:])/area # mask outside US, make units DALYs per km
plot_nsmokeDtotal = (D_nosmoke_total[1,:,:])/area
plot_Dtotal_smoke = (dD_total_smoke[1,:,:])/area
plot_afDsmoke = (dD_total_af[1,:,:])/area

# calc total US deaths for each so we can put on the map
US_Dtotal = np.nansum(D_total_total,axis=(1,2))
US_Dnosmoke = np.nansum(D_nosmoke_total,axis=(1,2))
US_Dsmoke = np.nansum(dD_total_smoke,axis=(1,2))
US_Dsmoke_af = np.nansum(dD_total_af,axis=(1,2))

total_str = str(round(US_Dtotal[0])) + '-' + str(round(US_Dtotal[-1]))
nosmoke_str = str(round(US_Dnosmoke[0])) + '-' + str(round(US_Dnosmoke[-1]))
smoke_str = str(round(US_Dsmoke[0])) + '-' + str(round(US_Dsmoke[-1]))
af_str = str(round(US_Dsmoke_af[0])) + '-' + str(round(US_Dsmoke_af[-1]))

fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(20, 13), constrained_layout=True,
                      subplot_kw={'projection': cart_proj})

axlist = axarr.flatten()
for ax in axlist:
    mk_map(ax)
c=axlist[0].pcolor(lon,lat,plot_Dtotal,shading='nearest',transform = ccrs.PlateCarree(),
        norm=colors.LogNorm(vmin=0.0001, vmax=10))
axlist[0].set_title('total PM DALYs '+total_str, fontsize=16)
cb1 = fig.colorbar(c, ax=axlist[0], orientation='horizontal', shrink=0.74, pad=0,
                   label = 'DALYs per km2 per year')
axlist[0].coastlines()

c=axlist[1].pcolor(lon,lat, plot_nsmokeDtotal, transform=ccrs.PlateCarree(),shading='nearest',
                norm=colors.LogNorm(vmin=0.0001, vmax=10))
axlist[1].set_title('total no smoke DALYs '+nosmoke_str, fontsize=16)
cb1 = fig.colorbar(c, ax=axlist[1], orientation='horizontal', shrink=0.74, pad=0,
                   label = 'DALYs per km2 per year')
axlist[1].coastlines()

c=axlist[2].pcolor(lon,lat, plot_Dtotal_smoke, transform=ccrs.PlateCarree(),shading='nearest',
                norm=colors.LogNorm(vmin=0.0001, vmax=1))
axlist[2].set_title('DALYs averted if no smoke '+smoke_str, fontsize=16)
cb = fig.colorbar(c, ax=axlist[2], orientation='horizontal', shrink=0.74, pad=0,
                  label = 'DALYs per km2 per year')

c=axlist[3].pcolor(lon,lat, plot_afDsmoke, transform=ccrs.PlateCarree(),shading='nearest',
                norm=colors.LogNorm(vmin=0.0001, vmax=1))
axlist[3].set_title('DALYs currently attributable to smoke '+af_str, fontsize=16)
cb = fig.colorbar(c, ax=axlist[3], orientation='horizontal', shrink=0.74, pad=0,
                  label = 'DALYs per km2 per year')

plt.savefig(out_fig_path + 'US_smoke_DALYs_'+desc+'.png')
plt.show()

#%% plot difference in non_smoke for both methods
'''
fig, axarr = plt.subplots(nrows=1, ncols=1, figsize=(20, 13), constrained_layout=True,
                      subplot_kw={'projection': ccrs.PlateCarree()})
## IMPORTANT - adjust lat lons to corners!!
plot_background(axarr)
diff = pm_nosmoke - pm_nosmoke_old
c=axarr.pcolormesh(pglon,pglat,diff, transform=ccrs.PlateCarree())
axarr.set_title('nonsmoke PM new-old', fontsize=16)
cb1 = fig.colorbar(c, ax=axarr, orientation='horizontal', shrink=0.74, pad=0,
                   label = 'nonsmoke PM new-old [ug m-3]')
plt.savefig(out_fig_path + 'compare_nosmoke_diff_'+desc+'.png')
plt.show()

## IMPORTANT - adjust lat lons to corners!!
fig, axarr = plt.subplots(nrows=1, ncols=1, figsize=(20, 13), constrained_layout=True,
                      subplot_kw={'projection': ccrs.PlateCarree()})
plot_background(axarr)
diff = pm_nosmoke - pm_nosmoke_old
pct_diff = 100.0*(diff/pm_nosmoke)
c=axarr.pcolormesh(pglon,pglat,pct_diff, transform=ccrs.PlateCarree())
axarr.set_title('nonsmoke PM new-old', fontsize=16)
cb1 = fig.colorbar(c, ax=axarr, orientation='horizontal', shrink=0.74, pad=0,
                   label = 'nonsmoke PM pct diff')
plt.savefig(out_fig_path + 'compare_nosmoke_pctdiff_'+desc+'.png')
plt.show()
'''