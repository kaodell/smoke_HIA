#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mk_state_grid.py
    python script to make a grid with assigned states to follow Kelsey's
    HIA code, a lot of state grid cells are added/modified by hand
Created on Thu Jun 18 09:28:29 2020
@author: kodell
"""
#%% user inputs
shapefile_fn_state = 'cb_2018_us_state_20m/cb_2018_us_state_20m'
shapefile_fn_country = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA scripts/Longitude_Graticules_and_World_Countries_Boundaries-shp/99bfd9e7-bb42-4728-87b5-07f8c8ac631c2020328-1-1vef4ev.lu5nk.shp'
pop_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/population/rg_population_2010.nc'
# csv with states and assoc. EPA regions ... ultimately not used in paper
epa_region_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/EPAregions.csv'

state_grid_out_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/state_grid'
EPAregion_grid_out_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/EPAregion_grid'
country_grid_out_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/country_grid'

#%% import modules
import numpy as np
import shapefile
import pylab as plb
import pandas as pd
import matplotlib as mplt
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

from wrf import (getvar, get_cartopy, cartopy_xlim,cartopy_ylim)
#%% load shapefiles needed and grid
pop_nc_fid = Dataset(pop_file)
glon = pop_nc_fid.variables['lon'][:].data
glat = pop_nc_fid.variables['lat'][:].data
pop = pop_nc_fid.variables['population'][:].data
pop_density = pop_nc_fid.variables['population_density'][:].data
pop_nc_fid.close()

#%% assign states to grid and EPA region
state_grid = np.empty(glon.shape).astype('str')
state_grid[:] = 'NA'
EPA_grid = np.empty(glon.shape).astype('int64')
EPA_grid[:] = 0

# load state shapefile and EPA regions grid
epa_regions = pd.read_csv(epa_region_fn)
epa_regions.set_index('State ID',inplace=True)

states_file = shapefile.Reader(shapefile_fn_state)
states_shp = states_file.shape(0)
state_records = states_file.records()
state_shapes = states_file.shapes()
si = 0
for j in range(len(state_records)):
        name = state_records[j][4]
        if name in ['HI','AK','PR']:
            continue
        print(name)
        plume_shp = state_shapes[j]
        for i in range(len(plume_shp.parts)):
            i0 = plume_shp.parts[i]
            if i < len(plume_shp.parts)-1:
            		i1 = plume_shp.parts[i+1] - 1
            else:
            		i1 = len(plume_shp.points)
            seg = plume_shp.points[i0:i1+1]
            mpath = mplt.path.Path(seg)
            points = np.array((glon.flatten(), glat.flatten())).T
            mask = mpath.contains_points(points).reshape(glon.shape)
            state_inds = np.where(mask==True)
            state_grid[state_inds] = name
            EPA_grid[state_inds] = epa_regions.loc[name]
            
        # add the coast
        if name in ['WA','CA','OR','FL']:
            coast_inds = np.array([state_inds[0],state_inds[1]-1])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]
            
            coast_inds = np.array([state_inds[0],state_inds[1]-2])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

        if name in ['GA','SC','NC','DE','MA','FL']:
            coast_inds = np.array([state_inds[0],state_inds[1]+1])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]
           
            coast_inds = np.array([state_inds[0],state_inds[1]+2])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]
        
        if name in ['RI']:
            coast_inds = np.array([state_inds[0]-1,state_inds[1]])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

            coast_inds = np.array([state_inds[0]-2,state_inds[1]])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

        if name in ['ME','AL','LA','FL']:
            coast_inds = np.array([state_inds[0]-1,state_inds[1]])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

            coast_inds = np.array([state_inds[0]-2,state_inds[1]])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

        MA_inds = ((139,138,135,134,134,135,135,134,136,138,140,140,140,140,134),
                   (303,305,303,306,300,304,305,305,305,304,301,302,303,304,301))
        state_grid[MA_inds[0],MA_inds[1]] = 'MA'                
        EPA_grid[MA_inds[0],MA_inds[1]] = epa_regions.loc['MA']

        VA_inds = ((103,102,104,97,98,97,105,99,100,101,102,102,103,103),
                   (281,281,279,283,282,284,276,283,283,283,283,281,281,285))
        state_grid[VA_inds[0],VA_inds[1]] = 'VA'                
        EPA_grid[VA_inds[0],VA_inds[1]] = epa_regions.loc['VA']


        if name in ['FL']:
            coast_inds = ((5,8,8,10,19,18),(263,274,275,277,265,265))
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

        if name in ['AL']:
            coast_inds = np.array([state_inds[0]-1,state_inds[1]])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]
        
        # fix alabama
        AL_inds = ((38,39,40),
                      (225,225,225))
        state_grid[AL_inds[0],AL_inds[1]] = 'AL'
        EPA_grid[AL_inds[0],AL_inds[1]] = epa_regions.loc['AL']
        
        # fix long island
        if name in ['NY']:
            coast_inds = ((140,139,138,137,137,137,136,133,130,129,128),
                          (270,269,268,265,262,261,259,259,258,257,256))
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

            # add long island
            state_grid[125:131,291:296] = name
            state_grid[124:128,287:290] = name
            state_grid[129,296] = name
            state_grid[130,297] = name
        
            EPA_grid[125:131,291:296] = epa_regions.loc[name]
            EPA_grid[124:128,287:290] = epa_regions.loc[name]
            EPA_grid[129,296] = epa_regions.loc[name]
            EPA_grid[130,297] = epa_regions.loc[name]

        # fix texas
        if name in ['TX']:
            coast_inds = ((17,26,27,29,29,45,46,44,43,42,23,25,15,5,30,31),
                          (172,183,184,178,186,189,189,190,191,191,179,182,171,172,188,190))               
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

        # fix CT
        if name in ['CT']:
            coast_inds = ((130,131,131),(292,295,296))
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

        # fix NJ
        if name in ['NJ']:
            coast_inds = ((117,118,119,120,121,122,115,116,112,113,114,123,123,123,113,113,112),
                        (289,289,289,289,289,289,288,288,287,287,287,289,288,287,285,284,285))
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

        # fix washington
        if name in ['WA']:
            WA_inds = ((171,174,172),(46,42,44))
            state_grid[WA_inds[0],WA_inds[1]] = name
            EPA_grid[WA_inds[0],WA_inds[1]] = epa_regions.loc[name]

        # fix great lakes region
        if name in ['MI']:
            coast_inds = np.array([state_inds[0]+1,state_inds[1]])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

            coast_inds = np.array([state_inds[0],state_inds[1]-1])
            state_grid[coast_inds[0],coast_inds[1]] = name
            EPA_grid[coast_inds[0],coast_inds[1]] = epa_regions.loc[name]

        MI_inds = ((135,136,143,148,148,148,147,146,153,153,137,145,144,153,142,119,118,122,
                    136,136,134,131,142,145,146,147,148,147,144,146,142,141,141,140,136,135,130,128,125,123,122),
                   (234,235,235,229,224,223,222,220,204,225,235,219,216,209,215,222,221,239,
                    237,239,240,241,235,232,230,229,228,227,226,223,224,224,223,222,221,221,221,222,223,223,223))
        state_grid[MI_inds[0],MI_inds[1]] = 'MI'
        EPA_grid[MI_inds[0],MI_inds[1]] = epa_regions.loc['MI']

        MN_inds = ((150,151,152,154,156,158,153,155),
                   (193,194,195,196,198,202,195,197))
        state_grid[MN_inds[0],MN_inds[1]] = 'MN'
        EPA_grid[MN_inds[0],MN_inds[1]] = epa_regions.loc['MN']

        # fix chicago
        state_grid[120,216] = 'IL'
        state_grid[118,217] = 'IL'
        
        EPA_grid[120,216] = epa_regions.loc['IL']
        EPA_grid[118,217] = epa_regions.loc['IL']

        if name in ['IN']:
            state_grid[117,218] = name
            state_grid[117,220] = name
            state_grid[117,219] = name
            
            EPA_grid[117,218] = epa_regions.loc[name]
            EPA_grid[117,220] = epa_regions.loc[name]
            EPA_grid[117,219] = epa_regions.loc[name]


        if name in ['OH']:
            OH_ex_inds = ((124,123,120,119,120,121,121,120),
                          (251,249,245,243,242,242,240,240))
            state_grid[OH_ex_inds] = name
            EPA_grid[OH_ex_inds] = epa_regions.loc[name]

        if name in ['PA']:
            state_grid[127,255] = name
            EPA_grid[127,255] = epa_regions.loc[name]

        # fix top part of FL panhandle for GA
        GA_inds = ((46,46,45,47,48),
                       (244,243,244,259,259))
        state_grid[GA_inds[0],GA_inds[1]] = 'GA'
        EPA_grid[GA_inds[0],GA_inds[1]] = epa_regions.loc['GA']

        # fix top part of FL panhandle for AL
        AL_inds = ((44,44,43,43,42,42,41,41,40,40),
                   (226,227,226,227,227,228,227,228,227,228))
        state_grid[AL_inds[0],AL_inds[1]] = 'AL'
        EPA_grid[AL_inds[0],AL_inds[1]] = epa_regions.loc['AL']
        
        # fix maryland
        MD_inds = ((105,108,108,109,110,111,112,107,108,105,104),
                   (281,285,278,278,278,278,278,279,279,280,283))
        state_grid[MD_inds[0],MD_inds[1]] = 'MD'
        EPA_grid[MD_inds[0],MD_inds[1]] = epa_regions.loc['MD']  
        
        # fix north carolina
        NC_inds = ((77,76,76,75,75,74,74,78,78,77),
                   (270,271,272,273,274,274,275,261,262,262))
        state_grid[NC_inds[0],NC_inds[1]] = 'NC'
        EPA_grid[NC_inds[0],NC_inds[1]] = epa_regions.loc['NC']  

        # wisconnsin
        WI_inds = ((131,130,125,133,140,139,138,137,136,139,138,141,142,143,141,141,
                    141,148,147,149,150,151,150,149),
                   (215,215,215,215,213,214,214,213,213,215,216,213,212,212,217,217,
                    216,200,200,198,198,197,195,193))
        state_grid[WI_inds[0],WI_inds[1]] = 'WI'
        EPA_grid[WI_inds[0],WI_inds[1]] = epa_regions.loc['WI'] 
        
        # lousiana
        state_grid[32,218] = 'LA'
        EPA_grid[32,218] = epa_regions.loc['LA'] 
        
        # mississippi
        state_grid[(38,38,37),(218,219,217)] = 'MS'
        EPA_grid[(38,38,37),(218,219,217)] = epa_regions.loc['MS']       

# save to files
np.save(state_grid_out_fn,state_grid)
np.save(EPAregion_grid_out_fn,EPA_grid)


#%% make country grid 
# don't really use this in the paper
# load world map shape file 
country_grid = np.empty(glon.shape).astype('str')
country_grid[:] = 'NA'
country_file = shapefile.Reader(shapefile_fn_country)
country_shp = country_file.shape(0)
country_records = country_file.records()
country_shapes = country_file.shapes()
si = 0

# lets assign all the added points for the states above to the US
US_inds = np.where(state_grid != 'NA')
country_grid[US_inds] = 'United States'

# save grid
np.save(country_grid_out_fn,country_grid)

print('data saved, making figures')

#%% plot to check
# state grid
# wrf plotting info: https://wrf-python.readthedocs.io/en/latest/plot.html 
# load wrf file to get cartopy projection
ncfile = Dataset('/Users/kodell/Local Google Drive /wrfout_d01_2016-07-22_00:00:00')
slp = getvar(ncfile,"slp")
cart_proj = get_cartopy(slp)
ncfile.close()

# Figure 1 - check specific state assignment
# plot data, pick a state, any state! (except AK or HI)
state_plot = 'MI'
state_plot_inds = np.where(state_grid==state_plot)
#make figure
fig = plt.figure()
ax = plt.axes(projection=cart_proj)
# Download and add the states and coastlines
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8)
plt.scatter(glon[state_plot_inds], glat[state_plot_inds],
            color='black',s=2,transform=crs.PlateCarree())
plt.title('state grid')   

# Figure 2 - plot population outside state to see if we're missing coasts
non_state_inds = np.where(state_grid!=state_plot)
fig = plt.figure()
ax = plt.axes(projection=cart_proj)
# Download and add the states and coastlines
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8)
plt.scatter(glon[non_state_inds], glat[non_state_inds],s=6,
           c=pop_density[non_state_inds],transform=crs.PlateCarree(),
           vmin=0,vmax=0.00001,cmap='viridis_r')
plt.scatter(glon[122,238], glat[122,238],
            color='red',s=4,transform=crs.PlateCarree())
# set bounds
plt.title('non-state pop')
ax.set_xlim(cartopy_xlim(slp[state_plot_inds[0].min()-2:state_plot_inds[0].max()+2,
                             state_plot_inds[1].min()-2:state_plot_inds[1].max()+2]))
ax.set_ylim(cartopy_ylim(slp[state_plot_inds[0].min()-2:state_plot_inds[0].max()+2,
                             state_plot_inds[1].min()-2:state_plot_inds[1].max()+2]))

# Figure 3 - plot population outside US to see if we're missing coasts
nonUS_inds = np.where(state_grid == 'NA')
plt.title('non-US pop')
fig = plt.figure()
ax = plt.axes(projection=cart_proj)
# Download and add the states and coastlines
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8)
plt.scatter(glon[nonUS_inds], glat[nonUS_inds],s=3,
           c=pop_density[nonUS_inds],transform=crs.PlateCarree(),
           vmin=0,vmax=0.00001,cmap='viridis_r')
