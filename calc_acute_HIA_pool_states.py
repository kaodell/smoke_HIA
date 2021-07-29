#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_acute_HIA_pool_state.py
    python script to calculate asthma ED by state
    
    v2 - 01.04.21 - updated to look at acute outcomes by EPA region and season
    01.15.20 - added custom grid option
    v3
    03.25.21 - altered to make annual maps by state, like for mortalities, but for each year and season
                copied to new code and renamed calc_aute_HIA_pool_states.py
Created on Wed Sep 30 16:48:03 2020
@author: kodell
"""
#%% user inputs
# file locations
acute_HIA_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/acute_outcomes/acute_events__byBreyregion_kPM_v2_medbk_final_final.npz'
state_grid_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/state_grid.npy'
pop_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/population/rg_population_2010.nc'
br_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/acute_baseline_rates_US.csv'

# out file and figure path and descrption
out_fp = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/acute_outcomes/'
out_fig_path ='/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/acute_outcomes/'
out_desc = '_byState_'+acute_HIA_file[-28:-4]+'_reviewer_edits'

# years
years = [2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
       2017, 2018]

# set figure fonts
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_family('sansserif')
font.set_name('Tahoma')
#%% import modules
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from wrf import (getvar, get_cartopy, cartopy_xlim,cartopy_ylim)

#%% load data
# acute HIA
np_fid = np.load(acute_HIA_file)
annual_event_all = np_fid['annual_event_all']
byseason_all = np_fid['byseason_all']
lons = np_fid['lons']
lats = np_fid['lats']
np_fid.close()

# state grid
state_grid = np.load(state_grid_file)

# population data
nc_fid = Dataset(pop_file, 'r')
population = nc_fid.variables['population'][:].data
nc_fid.close()

# baseline rates
br = pd.read_csv(br_file)
br.set_index('outcome',inplace=True)

#%% calc total of events by state and season
annual_event_sum = np.nansum(annual_event_all,axis=(2,3))

annual_event_sum_state = np.empty([annual_event_sum.shape[0],annual_event_sum.shape[1],49])
annual_event_sum_state[:] = 0.0

byseason_sum_state = np.empty([annual_event_sum.shape[0],annual_event_sum.shape[1],4,49])
byseason_sum_state[:] = 0.0

total_asth_ED_state = np.empty([1,annual_event_sum.shape[1],49])
total_asth_ED_state[:] = 0.0

states = np.unique(state_grid)
grid_asthED = np.empty([1,13,state_grid.shape[0],state_grid.shape[1]])
grid_asthED[:] = np.nan
grid_season_asthED = np.empty([1,13,4,state_grid.shape[0],state_grid.shape[1]])
grid_season_asthED[:] = np.nan

pct_asthED = np.empty([1,13,state_grid.shape[0],state_grid.shape[1]])
pct_asthED[:] = np.nan
pct_season_asthED = np.empty([1,13,4,state_grid.shape[0],state_grid.shape[1]])
pct_season_asthED[:] = np.nan

state_names = np.delete(states,np.where(states=='NA')[0])
for i in range(49):
    state_name = state_names[i]
    state_inds = np.where(state_grid==state_name)
    annual_event_sum_state[:,:,i] = np.nansum(annual_event_all[:,:,state_inds[0],
                                                            state_inds[1]],axis=(2))
    byseason_sum_state[:,:,:,i] = np.nansum(byseason_all[:,:,:,state_inds[0],
                                                            state_inds[1]],axis=(3))
    total_asth_ED_state[:,:,i] = np.nansum(population[state_inds[0],state_inds[1]])*br.loc['asth ED']['baseline rate']
    for j in range(13):
        # make gridded version for plotting
        grid_asthED[:,j,state_inds[0],state_inds[1]] = annual_event_sum_state[0,j,i]
        pct_asthED[:,j,state_inds[0],state_inds[1]] = 100.0*(annual_event_sum_state[0,j,i]/total_asth_ED_state[0,0,i])

        for k in range(4):
            grid_season_asthED[:,j,k,state_inds[0],state_inds[1]] = byseason_sum_state[0,j,k,i]
            pct_season_asthED[:,j,k,state_inds[0],state_inds[1]] = 100.0*(byseason_sum_state[0,j,k,i]/total_asth_ED_state[0,0,i])

#%% make figures

def mk_map(ax):
    ax.patch.set_visible(False)
    shapename = 'admin_1_states_provinces_lakes_shp'
    states_shp = shpreader.natural_earth(resolution='110m',category='cultural', name=shapename)
    facecolor = (0, 0, 0, 0)
    edgecolor = 'black'
    for state in shpreader.Reader(states_shp).geometries():
        ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)
# 1) events by state annually
'''
fig,axarr_g = plt.subplots(ncols=3,nrows=5,
                      subplot_kw={'projection': ccrs.PlateCarree()})
axarr=axarr_g.flatten()
for i in range(13):
    ax = axarr[i]
    mk_map(ax)
    ax.outline_patch.set_edgecolor('white')
    cs = ax.pcolor(lons,lats, pct_asthED[0,i,:,:],transform=ccrs.PlateCarree(),shading='nearest',cmap='Reds',
                   vmin=0,vmax=2.00,)
    ax.set_title(years[i],fontsize=10,fontproperties = font)
axarr[13].outline_patch.set_edgecolor('white')
axarr[14].outline_patch.set_edgecolor('white')
cax,kw = matplotlib.colorbar.make_axes(axarr[13],location='bottom',pad=-5)
cbar=fig.colorbar(cs,cax=cax,**kw,extend='max')
cbar.set_label('percent annual ED visits',fontsize=12,fontproperties = font)
plt.subplots_adjust(top=0.95,bottom=0.01,left=0.05,right=0.99,hspace=0.2,wspace=0.0)
plt.savefig(out_fig_path+'paper_annual_morbidity_map'+out_desc+'.png',dpi=500)
plt.show()
'''
# 2) by season 
si = 0
for season in ['JFM','AMJ','JAS','OND']:
    fig,axarr_g = plt.subplots(ncols=3,nrows=5,
                          subplot_kw={'projection': ccrs.PlateCarree()})
    axarr=axarr_g.flatten()
    for i in range(13):
        ax = axarr[i]
        mk_map(ax)
        ax.outline_patch.set_edgecolor('white')
        cs = ax.pcolor(lons,lats, pct_season_asthED[0,i,si,:,:],transform=ccrs.PlateCarree(),shading='nearest',cmap='Reds',
                       vmin=0,vmax=0.75,)
        ax.set_title(years[i],fontsize=10,fontproperties = font)
    axarr[13].outline_patch.set_edgecolor('white')
    axarr[14].outline_patch.set_edgecolor('white')
    cax,kw = matplotlib.colorbar.make_axes(axarr[13],location='bottom',pad=-5)
    cbar=fig.colorbar(cs,cax=cax,**kw,extend='max')
    cbar.set_label('% annual asthma ED visits in '+season,fontsize=12,fontproperties = font)
    plt.subplots_adjust(top=0.95,bottom=0.01,left=0.05,right=0.99,hspace=0.2,wspace=0.0)
    plt.savefig(out_fig_path+'paper_'+season+'_morbidity_map'+out_desc+'.png',dpi=500)
    plt.show()
    si += 1