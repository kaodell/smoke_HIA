#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
prep_smokePM_4HIA.py
    python script to load krigged PM from 2006-2018, calculate smokePM, 
    save output and make figures
Created on Wed Nov 11 14:55:17 2020
@author: kodell
"""
#%% user inputs
# file locations
pop_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/population/rg_population_2010.nc'
state_grid_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/state_grid.npy'
# check below for krigged PM files used
kpm_fp = '/Users/kodell/Desktop/data_not4gdrive/kriging/kriging data/kriging_v2/'
kpm_desc = '_v2_statfix_medbk'

# out file and figure path and descrption
out_file_path = '/Users/kodell/Desktop/data_not4gdrive/HIA_inputs/PM/'
out_fig_path ='/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/PMfigures/'
out_desc = '_v2_medbk_final'

# years to look at
years = [2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
       2017, 2018]

#%% import modules
import numpy as np
import plotly.graph_objects as go
from netCDF4 import Dataset
import netCDF4
import plotly.io as pio
from plotly.subplots import make_subplots
import pandas as pd
import datetime as dt
pio.renderers.default = "chrome"
pio.templates.default = "seaborn"

import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.feature import NaturalEarthFeature
from wrf import (getvar, get_cartopy, cartopy_xlim,cartopy_ylim)
    
#%% load data
# PM data
for year in years:
    # load kriging data
    # updated krige
    # can be downloaded at: https://doi.org/10.25675/10217/230602 
    kda_fid = Dataset(kpm_fp+'krigedPM25_'+str(year)+kpm_desc+'.nc')
    lon = kda_fid['lon'][:].data
    lat = kda_fid['lat'][:].data
    kda_PM1 = kda_fid['PM25'][:].data
    nosmokePM1 = kda_fid['Background_PM25'][:].data
    HMS_smoke = kda_fid['HMS_Smoke'][:].data
    kda_fid.close()

    print(year)

    # remove negatives from total and background; make zero
    # lowest background gets is -0.5; so not too negative; and it is rare.
    # but let's just remove it - doesn't make a difference overall
    nosmokePM = np.where(nosmokePM1<0,0,nosmokePM1)
    kda_PM = np.where(kda_PM1<0,0,kda_PM1)
    
    # test leaving negatives in ... doesn't impact results
    # kda_PM = kda_PM1
    # nosmokePM = nosmokePM1
    
    # calculate smoke PM
    smokePM = HMS_smoke*(kda_PM - nosmokePM)
        
    if year == years[0]:
        # create array
        # add nan for leap day so axis lines up
        ly_nans = np.empty([1,189,309])
        ly_nans[:] = np.nan
        smokePM_ly = np.insert(smokePM,59,ly_nans,axis=0)
        nosmokePM_ly = np.insert(nosmokePM,59,ly_nans,axis=0)
        totPM_ly = np.insert(kda_PM,59,ly_nans,axis=0)
        HMS_ly = np.insert(HMS_smoke,59,ly_nans,axis=0)

        allday_smokePM = np.array([smokePM_ly])
        allday_nosmokePM = np.array([nosmokePM_ly])
        allday_totPM = np.array([totPM_ly])
        allday_HMS = np.array([HMS_ly])

    elif year in [2008,2012,2016]:
        #these are leap years
        allday_smokePM = np.concatenate([allday_smokePM,np.array([smokePM])],axis=0)
        allday_nosmokePM = np.concatenate([allday_nosmokePM,np.array([nosmokePM])],axis=0)
        allday_totPM = np.concatenate([allday_totPM,np.array([kda_PM])],axis=0)
        allday_HMS = np.concatenate([allday_HMS,np.array([HMS_smoke])],axis=0)

    else:
        # add nan for leap day so axis lines up and stack
        smokePM_ly = np.insert(smokePM,59,ly_nans,axis=0)
        allday_smokePM = np.concatenate([allday_smokePM,np.array([smokePM_ly])],axis=0)
        
        nosmokePM_ly = np.insert(nosmokePM,59,ly_nans,axis=0)
        allday_nosmokePM = np.concatenate([allday_nosmokePM,np.array([nosmokePM_ly])],axis=0)
        
        totPM_ly = np.insert(kda_PM,59,ly_nans,axis=0)
        allday_totPM = np.concatenate([allday_totPM,np.array([totPM_ly])],axis=0)
        
        HMS_ly = np.insert(HMS_smoke,59,ly_nans,axis=0)
        allday_HMS = np.concatenate([allday_HMS,np.array([HMS_ly])],axis=0)

#%% save combined daily concentrations to netCDF file
outfn = out_file_path + 'krigedPM25_06-18' + out_desc + '.nc'
nc_w_fid = netCDF4.Dataset(outfn, 'w', format='NETCDF4',clobber=True)
nc_w_fid.description = 'Krigged PM2.5 concentrations from EPA AQS surface monitors for 2006-2018' + out_desc

# Need to define dimensions that will be used in the file
nc_w_fid.createDimension('year', len(years))
nc_w_fid.createDimension('DOY',allday_totPM.shape[1])
nc_w_fid.createDimension('lonx', lon.shape[0])
nc_w_fid.createDimension('lony', lon.shape[1])

# define variables for netcdf file
year_nc = nc_w_fid.createVariable('year', 'f8', ('year',))
lon_nc = nc_w_fid.createVariable('lon', 'f8', ('lonx','lony'))
lat_nc = nc_w_fid.createVariable('lat', 'f8', ('lonx','lony'))

PM25_nc = nc_w_fid.createVariable('PM25', 'f8', ('year', 'DOY', 'lonx','lony'))
PM25_background_nc = nc_w_fid.createVariable('Background PM25', 'f8', ('year', 'DOY', 'lonx', 'lony'))
PM25_smoke_nc = nc_w_fid.createVariable('Smoke PM25', 'f8', ('year', 'DOY', 'lonx', 'lony'))
HMS_smoke_nc = nc_w_fid.createVariable('HMS Smoke', 'f8', ('year', 'DOY', 'lonx', 'lony'))

year_nc.setncatts({'units':'date','long_name':'year',\
               'var_desc':'year'})
lat_nc.setncatts({'units':'degrees','long_name':'degrees latitude for data grid centers',\
               'var_desc':'Latitude (centers) [degrees]'})
lon_nc.setncatts({'units':'degrees','long_name':'degrees longitude for data grid centers',\
               'var_desc':'Longitude (centers) [degrees]'})

PM25_nc.setncatts({'units':'ug/m3','long_name':'PM2.5 concentration [ug/m3]',\
               'var_desc':'24 hour average PM2.5 concentration'})
PM25_background_nc.setncatts({'units':'ug/m3','long_name':'PM2.5 background concentration [ug/m3]',\
               'var_desc':'Seasonal PM25 Background (DJF, MAM, JJA, SON)'})
PM25_smoke_nc.setncatts({'units':'ug/m3','long_name':'Smoke-specific PM2.5',\
               'var_desc':'PM2.5 smoke concentration [ug/m3]'})
HMS_smoke_nc.setncatts({'units':'binary','long_name':'HMS Smoke',\
               'var_desc':'Binary HMS Smoke: 1 = smoke, 0 = no smoke'})

# data
year_nc[:] = years
lon_nc[:] = lon
lat_nc[:] = lat

PM25_nc[:] =  allday_totPM 
PM25_background_nc[:] =  allday_nosmokePM
PM25_smoke_nc[:] =  allday_smokePM
HMS_smoke_nc[:] = allday_HMS

nc_w_fid.close()
print('data saved, making figures')

#%% load population data
# load state grid to remove values outside the US
state_grid = np.load(state_grid_fn)
non_US_inds = np.where(state_grid=='NA')


nc_fid = Dataset(pop_file, 'r')
population = nc_fid.variables['population'][:].data
nc_fid.close()

#%% mask values outside the US and calc pop-weighted annual smoke PM
population[non_US_inds] = np.nan
allday_smokePM[:,:,non_US_inds[0],non_US_inds[1]]=np.nan
annavg_smokePM_ay = np.nanmean(allday_smokePM,axis=1)
pop_weighted_smokePM = np.nansum((annavg_smokePM_ay*population)/np.nansum(population),axis=(1,2))

allday_totPM[:,:,non_US_inds[0],non_US_inds[1]]=np.nan

# average smoke PM and total PM over full time period
cavg_smokePM = np.nanmean(allday_smokePM,axis=(0,1))
cavg_smokePM[non_US_inds] = np.nan

cavg_totPM = np.nanmean(allday_totPM,axis=(0,1))
cavg_totPM[non_US_inds] = np.nan

fraction = cavg_smokePM/cavg_totPM

#%% calc and print eUS and wUS pop-weighted smoke PM and population
pop_x_pm= annavg_smokePM_ay*population
# loop through states
east_pop = 0
west_pop = 0
west_popwPM = [0]*13
east_popwPM = [0]*13
eastPM = np.empty(annavg_smokePM_ay.shape)
eastPM[:] = np.nan
westPM = np.empty(annavg_smokePM_ay.shape)
westPM[:] = np.nan

for state in np.unique(state_grid):
    if state == 'NA':
        continue
    if state in ['WA','CA','OR','NV','AZ','NM','CO','UT','WY','MT','ID']:
        sinds = np.where(state_grid==state)
        west_pop += np.sum(population[sinds])
        west_popwPM += np.sum(pop_x_pm[:,sinds[0],sinds[1]],axis=(1))
        westPM = np.where(state_grid==state,annavg_smokePM_ay,westPM)
    else:
        sinds = np.where(state_grid==state)
        east_pop += np.sum(population[sinds])
        east_popwPM += np.sum(pop_x_pm[:,sinds[0],sinds[1]],axis=(1))
        eastPM = np.where(state_grid==state,annavg_smokePM_ay,eastPM)

print('wUS pop',west_pop,'eUS pop',east_pop)
print('pop-weighted smoke PM wUS',west_popwPM/west_pop,'eUS',east_popwPM/east_pop)

#%% make PM figures
# 1) populeation-weighted PM by year, Figure 1b
fig = go.Figure()
meanUSsmoke = np.nanmean(allday_smokePM,axis=(1,2,3))
fig.add_trace(go.Scatter(y=meanUSsmoke,x=years,name='mean US smoke PM<sub>2.5</sub>',
                         mode='markers',marker_size=5,marker_color='black'))
'''
fig.add_trace(go.Scatter(y=pop_weighted_smokePM,x=years,name='pop-weighted smoke PM<sub>2.5</sub>',
                         mode='markers',marker_size=5,marker_color='darkorange'))
'''
fig.add_trace(go.Scatter(y=east_popwPM/east_pop,x=years,name='eUS pop-weighted smoke PM<sub>2.5</sub>',
                         mode='markers',marker_size=5,marker_color='darkorange'))
fig.add_trace(go.Scatter(y=west_popwPM/west_pop,x=years,name='wUS pop-weighted smoke PM<sub>2.5</sub>',
                         mode='markers',marker_size=5,marker_color='red'))

# dist of avg smoke in each grid cell
gc_annavg_smoke = np.nanmean(allday_smokePM,axis=(1))
yi = 0
for year in years: 
    smoke_plot = gc_annavg_smoke[yi,:,:].flatten()
    x = [years[yi]]*len(smoke_plot)
    fig.add_trace(go.Box(y=smoke_plot,
                         x=x,showlegend=False,
                         marker_color='gray',boxpoints=False))
    yi += 1
fig.update_layout(yaxis_type='log',yaxis_range=(-2,1.5),plot_bgcolor='white',
                  title_text='b) Annual Mean Smoke PM<sub>2.5</sub>')
fig.update_xaxes(title_text='year',
                 ticks="outside", tickvals= [2006,2008,2010,2012,2014,2016,2018],
                 tickwidth=1, tickcolor='black', ticklen=10,
                 showline=True, linewidth=1, linecolor='black')
fig.update_yaxes(title_text='smoke PM<sub>2.5</sub> [&#956;gm<sup>-3</sup>]',
                 tickvals= [0.01,0.1,1.0,10.0],
                 ticks="outside", tickwidth=1, tickcolor='black', ticklen=10,
                 showline=True, linewidth=1, linecolor='black')
fig.write_image(out_fig_path + 'annual_PM_distribution'+out_desc+'.png',
                height=400,width=700,scale=4)
fig.show()

#%% maps of PM
# map plotting details
# wrf plotting info: https://wrf-python.readthedocs.io/en/latest/plot.html 
# load wrf file to get cartopy projection
ncfile = Dataset('/Users/kodell/Desktop/data_not4gdrive/wrfout_d01_2016-07-22_00_00_00')
slp = getvar(ncfile,"slp")
cart_proj = get_cartopy(slp)
ncfile.close()
del slp

def mk_map():
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0, 0, 1, 1], projection=cart_proj)
    plt.gca().outline_patch.set_visible(False)
    ax.patch.set_visible(False)
    shapename = 'admin_1_states_provinces_lakes_shp'
    states_shp = shpreader.natural_earth(resolution='110m',category='cultural', name=shapename)
    facecolor = (0, 0, 0, 0)
    edgecolor = 'black'
    for state in shpreader.Reader(states_shp).geometries():
        ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)
    return fig,ax


# 2) avg smoke PM over full time period, Figure 1a
fig, ax = mk_map()
cmap = matplotlib.cm.get_cmap('OrRd',10)
bounds = np.round(np.logspace(np.log10(0.03),np.log10(2.0),9,base=10.0),2)
norm=matplotlib.colors.BoundaryNorm(boundaries=bounds,ncolors=10,extend='max')
cs = ax.pcolor(lon,lat, cavg_smokePM,transform=ccrs.PlateCarree(),norm=norm,cmap=cmap,shading='nearest')
cbar=plt.colorbar(cs,location='bottom',shrink=0.60,pad=0.001)
cbar.set_label(r'Smoke PM$_{2.5}$ [$\mu$g m$^{-3}$]',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.title('a) Mean Smoke PM$_{2.5}$ 2006-2018',fontsize=20)
plt.savefig(out_fig_path+'avg_smokePM_06-18_reds_log'+out_desc+'.png',dpi=300)
plt.show()

# 3) map of total PM, Figure S2a
fig, ax = mk_map()
cmap = matplotlib.cm.get_cmap('OrRd',9)
bounds = [2,4,6,8,10,12,14,16,18]
norm=matplotlib.colors.BoundaryNorm(boundaries=bounds,ncolors=9)
cs = ax.pcolor(lon,lat, cavg_totPM,transform=ccrs.PlateCarree(),norm=norm,cmap=cmap,shading='nearest')
cbar=plt.colorbar(cs,location='bottom',shrink=0.60,pad=0.001)
cbar.set_label(r'PM$_{2.5}$ [$\mu$g m$^{-3}$]',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.title('a) Mean PM$_{2.5}$ 2006-2018',fontsize=18)
plt.savefig(out_fig_path+'avgPM_06-18_reds_log'+out_desc+'.png',dpi=300)
plt.show()

# 3) smoke fraction, Figure S2b
fig, ax = mk_map()
cmap = matplotlib.cm.get_cmap('Greys',6)
bounds = [0,0.01,0.025,0.05,0.1,0.25,1]
norm=matplotlib.colors.BoundaryNorm(boundaries=bounds,ncolors=6)
cs = ax.pcolor(lon,lat, fraction,transform=ccrs.PlateCarree(),norm=norm,cmap=cmap,shading='nearest')
cbar=plt.colorbar(cs,location='bottom',shrink=0.60,pad=0.001)
cbar.set_label(r'Smoke Fraction PM$_{2.5}$',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.title('b) Smoke Fraction 2006-2018',fontsize=18)
plt.savefig(out_fig_path+'smoke_fraction_PM_06-18_log'+out_desc+'.png',dpi=300)
plt.show()


#%% figure graveyard
'''
#%% additional figures
# make population figure
# calc population density
population_density = population.data/225.0

# set up a map
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
plt.gca().outline_patch.set_visible(False)
ax.patch.set_visible(False)
shapename = 'admin_1_states_provinces_lakes_shp'
states_shp = shpreader.natural_earth(resolution='110m',
                                         category='cultural', name=shapename)

for state in shpreader.Reader(states_shp).geometries():
    # pick a default color for the land with a black outline,
    # this will change if the storm intersects with our track
    facecolor = (0, 0, 0, 0)
    edgecolor = 'black'
    ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)

cs = ax.pcolor(glon.data,glat.data, population_density,transform=ccrs.PlateCarree(),shading='nearest',
           norm=matplotlib.colors.LogNorm(),vmin=1,vmax=1000,cmap='Blues')
cbar=plt.colorbar(cs,location='bottom',shrink=0.70,pad=0.01,extend='max')
cbar.set_label(r'people per km$^2$')
plt.title('US Population Density, 2015',fontsize=16)
plt.savefig(out_fig_path+'population_map.png',dpi=300)
plt.show()

#%% make pop x PM figure to check eUS pop-weighted somke
yi = 0
for year in years:
    year_avg = np.nanmean(allday_smokePM[yi,:,:,:],axis=(0))
    pop_PM = population*year_avg
    # set up a map
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
    plt.gca().outline_patch.set_visible(False)
    ax.patch.set_visible(False)
    shapename = 'admin_1_states_provinces_lakes_shp'
    states_shp = shpreader.natural_earth(resolution='110m',
                                             category='cultural', name=shapename)
    
    for state in shpreader.Reader(states_shp).geometries():
        # pick a default color for the land with a black outline,
        # this will change if the storm intersects with our track
        facecolor = (0, 0, 0, 0)
        edgecolor = 'black'
        ax.add_geometries([state], ccrs.PlateCarree(),
                          facecolor=facecolor, edgecolor=edgecolor)
    
    cs = ax.pcolor(glon.data,glat.data, pop_PM,transform=ccrs.PlateCarree(),shading='nearest',
                   norm=matplotlib.colors.LogNorm(),vmin=1.0)
    cbar=plt.colorbar(cs,location='bottom',shrink=0.70,pad=0.01,extend='max')
    cbar.set_label(r'Smoke PM$_{2.5}$ [$\mu$gm$^{-3}$] x population')
    plt.title('pop x smokePM ' +str(year),fontsize=16)
    plt.savefig(out_fig_path+'/PMfigures/pop_x_PM'+str(year)+'.png',dpi=300)
    #plt.show()
    yi += 1

#%% investigate eUS pm
NYinds = np.where(state_grid=='NY')
fig = go.Figure()
yi = -1
#for gi in range(len(NYinds[0])):
fig.add_trace(go.Scatter(y=np.nanmedian(allday_smokePM[yi,:,NYinds[0],NYinds[1]],axis=0),name='smoke PM2.5',
                         mode='lines',marker_size=0.5,marker_color='orange',
                         showlegend=False))
fig.add_trace(go.Scatter(y=np.nanmedian(allday_nosmokePM[yi,:,NYinds[0],NYinds[1]],axis=0),name='smoke PM2.5',
                         mode='lines',marker_size=0.5,marker_color='blue',
                         showlegend=False))

fig.show()

# HMS days map
HMS_sum = np.nansum(allday_HMS,axis=(0,1))
HMS_percent = 100.0*(HMS_sum/4755.0)

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
plt.gca().outline_patch.set_visible(False)
ax.patch.set_visible(False)
shapename = 'admin_1_states_provinces_lakes_shp'
states_shp = shpreader.natural_earth(resolution='110m',
                                         category='cultural', name=shapename)
facecolor = (0, 0, 0, 0)
edgecolor = 'grey'
for state in shpreader.Reader(states_shp).geometries():
    ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)
cmap = matplotlib.cm.get_cmap('magma',10)
bounds = (0,2,4,6,8,10,12,14,16,18)
norm=matplotlib.colors.BoundaryNorm(boundaries=bounds,ncolors=10)

cs = ax.pcolor(pglon.data,pglat.data, HMS_percent,transform=ccrs.PlateCarree(),
               norm=norm,cmap=cmap)
cbar=plt.colorbar(cs,location='bottom',shrink=0.60,pad=0.001)
cbar.set_label(r'percent smoke days')
plt.title('Percent Smoke Days 2006-2018',fontsize=14)
plt.savefig(out_fig_path+'pct_HMS_06-18'+out_desc+'.png',dpi=300)
plt.show()

# for checking PM details
# residual maps
residual = allday_totPM - allday_nosmokePM
residual_nonsmokedays = np.where(allday_HMS==0,residual,np.nan)

smokePM_smokedays = np.where(allday_HMS==1,allday_smokePM,np.nan)

# residual on non-smoke days
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
plt.gca().outline_patch.set_visible(False)
ax.patch.set_visible(False)
shapename = 'admin_1_states_provinces_lakes_shp'
states_shp = shpreader.natural_earth(resolution='110m',
                                         category='cultural', name=shapename)
facecolor = (0, 0, 0, 0)
edgecolor = 'grey'
for state in shpreader.Reader(states_shp).geometries():
    ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)
cmap = matplotlib.cm.get_cmap('OrRd',10)
bounds = np.round(np.logspace(np.log10(0.03),np.log10(2.0),9,base=10.0),2)
norm=matplotlib.colors.BoundaryNorm(boundaries=bounds,extend='max',ncolors=10)

cs = ax.pcolor(pglon.data,pglat.data, np.nanmean(residual_nonsmokedays,axis=(0,1)),
               transform=ccrs.PlateCarree(),norm=norm,cmap=cmap)
cbar=plt.colorbar(cs,location='bottom',shrink=0.60,pad=0.001)
cbar.set_label(r'residual PM2.5')
plt.title('Residual PM2.5 on nonsmoke days 2006-2018',fontsize=14)
plt.savefig(out_fig_path+'residual_06-18'+out_desc+'.png',dpi=300)
plt.show()


# smoke PM on smoke days
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
plt.gca().outline_patch.set_visible(False)
ax.patch.set_visible(False)
shapename = 'admin_1_states_provinces_lakes_shp'
states_shp = shpreader.natural_earth(resolution='110m',
                                         category='cultural', name=shapename)
facecolor = (0, 0, 0, 0)
edgecolor = 'grey'
for state in shpreader.Reader(states_shp).geometries():
    ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)
cmap = matplotlib.cm.get_cmap('OrRd',10)
bounds = np.round(np.logspace(np.log10(1.0),np.log10(10),9,base=10.0),2)
norm=matplotlib.colors.BoundaryNorm(boundaries=bounds,extend='max',ncolors=10)

cs = ax.pcolor(pglon.data,pglat.data, np.nanmean(smokePM_smokedays,axis=(0,1)),
               transform=ccrs.PlateCarree(),norm=norm,cmap=cmap)
cbar=plt.colorbar(cs,location='bottom',shrink=0.60,pad=0.001)
cbar.set_label(r'smoke PM2.5 on smoke days')
plt.title('Smoke PM2.5 on smoke days 2006-2018',fontsize=14)
plt.savefig(out_fig_path+'smoke_on_smokedays_06-18'+out_desc+'.png',dpi=300)
plt.show()
'''
