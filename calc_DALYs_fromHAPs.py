#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
calc_DALYs_fromHAPs.py

a python script to load PM smoke exposure, multiply by
HAPs to PM ratios, then use equation to estimate DALYs

this code has now been moved to a new project folder, smoke-specific HIA
and has been expanded to work with all HAPs with a PM relationship in WE-CAN

Created on Thu Nov 14 11:35:19 2019
@author: kodell
"""
#%% Load Modules
import pandas as pd
import numpy as np
import scipy.stats as st
from netCDF4 import Dataset

import plotly
import plotly.graph_objs as go
import plotly.io as pio
pio.renderers.default = "chrome"
pio.templates.default = "seaborn"

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import cartopy.io.shapereader as shpreader
from wrf import (getvar, get_cartopy, cartopy_xlim,cartopy_ylim)
#%% User Inputs
years = np.arange(2006,2019)

# kriging data and HMS data paths
kPM_fp = '/Users/kodell/Desktop/data_not4gdrive/HIA_inputs/PM/'
kPM_desc = '_v2_medbk_final' # which kPM file to load
PM_DALY_desc = kPM_desc+'_pop2010'

# ratio file name
# these ratios are in ug m-3 HAP/ ug m-3 PM and can be downloaded
# from the supplementary material of O'Dell et al. (2020) https://doi.org/10.1021/acs.est.0c04497
ratio_data_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/VOC2PM_ratios_4paper_NASAmrg_R1TOGAupdate_nonpar_stats_2.5_97.5_test.csv'

# DALY damage and effect factors file
# these are the combined damage and effect factors from Huijbregts et al. (2005) 
# appendix 1, https://doi.org/10.1897/2004-007R.1
DALY_EF_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/HAPs/DALY_effect_factors.csv'

# population file
pop_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/population/rg_population_2010.nc'

# state grid assignment file
state_grid_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/state_grid.npy'

# PM DALYs file
pmDALY_fn ='/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/mortalities/wilfire_smoke_PM_DALYs_gemm'+PM_DALY_desc+'.nc'

# out figure path and description
fig_path =  '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/DALYs/'
fig_desc = kPM_desc+'_2010pop'

#%% user-defined functions
def make_map(cart_projection):
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection=cart_proj)
    plt.gca().outline_patch.set_visible(False)
    ax.patch.set_visible(False)
    shapename = 'admin_1_states_provinces_lakes_shp'
    states_shp = shpreader.natural_earth(resolution='110m',category='cultural', name=shapename)
    facecolor = (0, 0, 0, 0)
    edgecolor = 'black'
    for state in shpreader.Reader(states_shp).geometries():
        ax.add_geometries([state], crs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)
    return fig, ax

#%% Load Files
# 1) Load kPM 
# PM data
kPM_fid = Dataset(kPM_fp + 'krigedPM25_06-18'+kPM_desc+'.nc')
smokePM_daily_allyears = kPM_fid['Smoke PM25'][:].data
glon = kPM_fid['lon'][:].data
glat = kPM_fid['lat'][:].data
kPM_fid.close()

# 2) Load TOGA mrg VOC:PM relationships
VOC2PMstats = pd.read_csv(ratio_data_fn)

# 3) load population
nc_fid = Dataset(pop_file, 'r')
population = nc_fid.variables['population'][:]
nc_fid.close()

# 4) load DALY/kg HAP numbers
DALY_effect_factors = pd.read_csv(DALY_EF_fn)
HAPs = DALY_effect_factors['Substance'].values

# 5) load PM DALYs
nc_fid = Dataset(pmDALY_fn)
PM_DALYs = nc_fid['smoke_DALY_fraction'][:]
tot_PM_DALYs = nc_fid['total_DALY'][:]
nc_fid.close()

#%% calc 06-18 avg smoke PM over US
state_grid = np.load(state_grid_fn)

# now, calculate average smokePM 2006-2018
decadal_avg_smokePM = np.nanmean(smokePM_daily_allyears,axis=(0,1))

nonUSinds = np.where(state_grid =='NA')
decadal_avg_smokePM[nonUSinds] = np.nan

#%% prep data for smoke

# let's use the young-age for all parts of the US (but test different ones)
VOC2PMstats.set_index('VOC name',inplace=True)
ages = ['young VOC PM1 ratio, median [ugm-3 STP/ugm-3 STP]',
        'medium VOC PM1 ratio, median [ugm-3 STP/ugm-3 STP]',
        'old VOC PM1 ratio, median [ugm-3 STP/ugm-3 STP]']

DALY_effect_factors.set_index('Substance',inplace=True)

# volume air inhaled per year 
V_air_yr = 14.4*365.0 #m3/yr ... average daily air intake from Logue et al. (2012)

# 95% x values are within a factor of k
DALY_effect_factors['carcinogenic DALY/kg uci']=DALY_effect_factors['carcinogenic DALY/kg'].values * DALY_effect_factors['carcinogeric k∂D/∂I'].values
DALY_effect_factors['carcinogenic DALY/kg lci']=DALY_effect_factors['carcinogenic DALY/kg'].values / DALY_effect_factors['carcinogeric k∂D/∂I'].values

DALY_effect_factors['non-carcinogenic DALY/kg uci']=DALY_effect_factors['non-carcinogenic DALY/kg'].values * DALY_effect_factors['non-carcinogeric k∂D/∂I'].values
DALY_effect_factors['non-carcinogenic DALY/kg lci']=DALY_effect_factors['non-carcinogenic DALY/kg'].values / DALY_effect_factors['non-carcinogeric k∂D/∂I'].values

#%% calc DALYs
dec_DALY = np.empty([3,len(HAPs),189,309])
dec_DALY[:] = np.nan

dec_DALY_uci = np.empty([3,len(HAPs),189,309])
dec_DALY_uci[:] = np.nan

dec_DALY_lci = np.empty([3,len(HAPs),189,309])
dec_DALY_lci[:] = np.nan

j = 0
for age in ages:
    i = 0
    for HAP in HAPs:
        # calc DALY concentration in ug/m3 and convert to kg/m3
        dec_Ckg = VOC2PMstats.loc[HAP][age]*decadal_avg_smokePM*(10.0**-9)
      
        # pull DALYs per kg of the HAP
        dDdI_c = float(DALY_effect_factors.loc[HAP]['carcinogenic DALY/kg'])
        dDdI_nc = float(DALY_effect_factors.loc[HAP]['non-carcinogenic DALY/kg'])
        
        dDdI_c_uci = float(DALY_effect_factors.loc[HAP]['carcinogenic DALY/kg uci'])
        dDdI_nc_uci = float(DALY_effect_factors.loc[HAP]['non-carcinogenic DALY/kg uci'])

        dDdI_c_lci = float(DALY_effect_factors.loc[HAP]['carcinogenic DALY/kg lci'])
        dDdI_nc_lci = float(DALY_effect_factors.loc[HAP]['non-carcinogenic DALY/kg lci'])

        # calc DALYs
        dec_DALY_c = dec_Ckg*V_air_yr*dDdI_c*population
        dec_DALY_nc = dec_Ckg*V_air_yr*dDdI_nc*population
        dec_DALY[j,i,:,:] = np.nansum([dec_DALY_c, dec_DALY_nc],axis=0)
        
        dec_DALY_c_uci = dec_Ckg*V_air_yr*dDdI_c_uci*population
        dec_DALY_nc_uci = dec_Ckg*V_air_yr*dDdI_nc_uci*population
        dec_DALY_uci[j,i,:,:] = np.nansum([dec_DALY_c_uci, dec_DALY_nc_uci],axis=0)

        dec_DALY_c_lci = dec_Ckg*V_air_yr*dDdI_c_lci*population
        dec_DALY_nc_lci = dec_Ckg*V_air_yr*dDdI_nc_lci*population
        dec_DALY_lci[j,i,:,:] = np.nansum([dec_DALY_c_lci, dec_DALY_nc_lci],axis=0)
        i+=1
    j+=1

dec_DALY_all = np.nansum(dec_DALY[0,:,:,:],axis=0)
dec_DALY_all_uci = np.nansum(dec_DALY_uci[0,:,:,:],axis=0)
dec_DALY_all_lci = np.nansum(dec_DALY_lci[0,:,:,:],axis=0)

#%% make plots
#  # wrf plotting info: https://wrf-python.readthedocs.io/en/latest/plot.html 
# load wrf file to get cartopy projection
ncfile = Dataset('/Users/kodell/Local Google Drive /wrfout_d01_2016-07-22_00:00:00')
slp = getvar(ncfile,"slp")
cart_proj = get_cartopy(slp)
ncfile.close()

# 1) average acrolein exposure
USdec_acrolein_Cug = VOC2PMstats.loc['acrolein'][ages[0]]*decadal_avg_smokePM
fig, ax = make_map(cart_proj)
cmap = matplotlib.cm.get_cmap('Blues',9)
bounds = np.round(np.logspace(np.log10(0.0004),np.log10(0.02),9,base=10.0),4)
norm=matplotlib.colors.BoundaryNorm(boundaries=bounds,ncolors=9,extend='max')
cs = ax.pcolor(glon, glat, USdec_acrolein_Cug,transform=crs.PlateCarree(),
                shading = 'nearest',norm=norm,cmap=cmap)
cbar=plt.colorbar(cs,location='bottom',shrink=0.60,pad=0.001)
cbar.set_label(r'Acrolein [$\mu$gm$^{-3}$]')
plt.title('13-yr Average Acrolein Exposure from Wildfire Smoke')
plt.savefig(fig_path + 'Acrolein_conc'+fig_desc+'_fullUS.png',dpi=300)
plt.show()

# 2) all DALYs from HAPs map
fig, ax = make_map(cart_proj)
cs=plt.pcolor(glon, glat, dec_DALY_all,transform=crs.PlateCarree(),
                 vmin=0,vmax=.01,shading = 'nearest') # take grid centers
plt.title('DALYs for 13-yr HAPs in Wildfire Smoke')
plt.colorbar(cs,extend='max',label = r'DALYs per year')
plt.savefig(fig_path + 'DALY_map'+fig_desc+'_fullUS.png',dpi=300)
plt.show()

# 3) US total DALYs from indv HAPs and PM
fig = go.Figure()
# sort HAPs by decreasing young DALYs
tot_DALYs = np.nansum(dec_DALY,axis=(2,3))  
tot_DALYs_uci = np.nansum(dec_DALY_uci,axis=(2,3))                         
tot_DALYs_lci = np.nansum(dec_DALY_lci,axis=(2,3))                         
                       
tot_DALYs_df = pd.DataFrame(data={'VOC':HAPs,'young':tot_DALYs[0,:],
                                  'young lci':tot_DALYs_lci[0,:],'young uci':tot_DALYs_uci[0,:]})                         
sorted_tot_DALYs_df = tot_DALYs_df.sort_values(by='young',ascending=True)


# note [14:] is to remove HAPs with no values on this figure (for supplement figure)
# note [22:] is to remove HAPs with DALYs < 0.01 (for main text figure)

fig.add_trace(go.Scatter(name='smoke HAPs', 
                         y=sorted_tot_DALYs_df['VOC'].values, 
                         x=sorted_tot_DALYs_df['young'].values,
                         mode='markers',
                         text = sorted_tot_DALYs_df['VOC'].values,
                         error_x = dict(type='data',symmetric=False,
                                        array = sorted_tot_DALYs_df['young uci'].values - sorted_tot_DALYs_df['young'].values,
                                        arrayminus = sorted_tot_DALYs_df['young'].values - sorted_tot_DALYs_df['young lci'].values),
                         textposition='middle right',
                         marker_size=10,marker_color='dimgray',showlegend=False)),

fig.add_trace(go.Scatter(name='smoke HAPs', y=['all HAPs'], x=[np.nansum(dec_DALY[0,:,:,:])],
                         mode='markers',text=['all HAPs'],textposition='middle right',
                         error_x = dict(type='data',symmetric=False,
                                        array = [np.nansum(dec_DALY_uci[0,:,:,:]) - np.nansum(dec_DALY[0,:,:,:])],
                                        arrayminus = [np.nansum(dec_DALY[0,:,:,:]) - np.nansum(dec_DALY_lci[0,:,:,:])]),
                         marker_size=10,marker_color='black',showlegend=False)),

fig.add_trace(go.Scatter(name='smoke PM2.5', y=['PM<sub>2.5</sub>'], x=[np.nansum(PM_DALYs[1,:,:])],
                         mode='markers',text=['smoke PM2.5'],textposition='middle right',
                        error_x=dict(type='data',symmetric=False,
                                      array=[np.nansum(PM_DALYs[2,:,:])-np.nansum(PM_DALYs[1,:,:])],
                                      arrayminus=[np.nansum(PM_DALYs[1,:,:])-np.nansum(PM_DALYs[0,:,:])]),
                        marker_size=10,marker_color='black',marker_symbol='square',showlegend=False))

fig.update_xaxes(type='log',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 showline=True, linewidth=1, linecolor='black',title='DALYs y<sup>-1</sup>')
fig.update_yaxes(ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 showline=True, linewidth=1, linecolor='black')
fig.update_layout(title='DALYs from HAPs and PM in US wildfire smoke, 2006-2018',
                  plot_bgcolor='rgba(0,0,0,0)',legend_traceorder='reversed',
                  font_family='Arial')
fig.write_image('/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/DALYs/PMandHAP0618DALYs_'+fig_desc+'_rotaxes_HAPs_all.png',
                scale=4,height=600,width=600)
fig.show()

#%% old plot graveyard
'''
i = 1
for year in years:
    fig = pl.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # setup of basemap ('lcc' = lambert conformal conic).
    # use major and minor sphere radii from WGS84 ellipsoid.
    m = Basemap(llcrnrlon=-125.5,llcrnrlat=25.,urcrnrlon=-65,urcrnrlat=50.,
                resolution='l',area_thresh=1000.,projection='merc',lon_0=-107.)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    
    # plot data
    cs = m.pcolormesh(pglon,pglat, anavg_smokePM_allyears[i,:,:],latlon=True,vmax=2)
    # add colorbar
    m.colorbar(cs,"bottom", size="5%", pad='2%')
    ax.set_title('Annual Average Smoke PM' + str(years[i-1]))
    pl.show()
    i += 1


fig = go.Figure()
fig.add_trace(go.Scatter(name='smoke PM2.5', x=['PM'], y=[np.nansum(wUS_PM_DALYs[1,:,:])],
                         error_y=dict(type='data',symmetric=False,
                                      array=[np.nansum(wUS_PM_DALYs[2,:,:])-np.nansum(wUS_PM_DALYs[1,:,:])],
                                      arrayminus=[np.nansum(wUS_PM_DALYs[1,:,:])-np.nansum(wUS_PM_DALYs[0,:,:])]),
                         marker_size=10,marker_color='steelblue'))
                         
fig.add_trace(go.Scatter(name='young smoke HAPs', x=['all HAPs'], y=[np.nansum(wUSdec_DALY[0,:,:,:])],
                         marker_size=10,marker_color='purple')),
fig.add_trace(go.Scatter(name='medium smoke HAPs', x=['all HAPs'], y=[np.nansum(wUSdec_DALY[1,:,:,:])],
                         marker_size=10,marker_color='red')),
fig.add_trace(go.Scatter(name='old smoke HAPs', x=['all HAPs'], y=[np.nansum(wUSdec_DALY[2,:,:,:])],
                         marker_size=10,marker_color='orange')),

# sort HAPs by decreasing young DALYs
wUStot_DALYs = np.nansum(wUSdec_DALY,axis=(2,3))                         
wUStot_DALYs_df = pd.DataFrame(data={'VOC':HAPs,'young':wUStot_DALYs[0,:],
                                  'medium':wUStot_DALYs[1,:],'old':wUStot_DALYs[2,:]})                         
wUSsorted_tot_DALYs_df = wUStot_DALYs_df.sort_values(by='young',ascending=False)
    
fig.add_trace(go.Scatter(name='young smoke HAPs', 
                         x=wUSsorted_tot_DALYs_df['VOC'].values, 
                         y=wUSsorted_tot_DALYs_df['young'].values,mode='markers',
                         marker_size=10,marker_color='purple',showlegend=False)),
fig.add_trace(go.Scatter(name='medium smoke HAPs', 
                         x=wUSsorted_tot_DALYs_df['VOC'].values, 
                         y=wUSsorted_tot_DALYs_df['medium'].values,mode='markers',
                         marker_size=10,marker_color='red',showlegend=False)),
fig.add_trace(go.Scatter(name='old smoke HAPs', 
                         x=wUSsorted_tot_DALYs_df['VOC'].values, 
                         y=wUSsorted_tot_DALYs_df['old'].values,mode='markers',
                         marker_size=10,marker_color='orange',showlegend=False)),

fig.update_yaxes(type='log')
fig.update_layout(title='DALYs from HAPs and PM in western US wildfire smoke, 2006-2018',
                  plot_bgcolor='rgba(0,0,0,0)')
fig.write_image('/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/DALYs/PMandHAP0618DALYs_wUS_'+fig_desc+'.png',
                scale=4,height=600,width=900)
fig.show()


i = 0
for year in years:
    fig = pl.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # setup of basemap ('lcc' = lambert conformal conic).
    # use major and minor sphere radii from WGS84 ellipsoid.
    m = Basemap(llcrnrlon=-125.5,llcrnrlat=25.,urcrnrlon=-65,urcrnrlat=50.,
            resolution='l',area_thresh=1000.,projection='merc',lon_0=-107.)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    
    # plot data
    cs = m.pcolormesh(pglon,pglat, annAcrolein[i,:,:],latlon=True,vmin=0)
    # add colorbar
    m.colorbar(cs,"bottom", size="5%", pad='2%',label = 'Acrolein [ugm-3]')
    ax.set_title('Acrolein from Smoke Exposure '+str(year))
    pl.show()
    i += 1
'''


