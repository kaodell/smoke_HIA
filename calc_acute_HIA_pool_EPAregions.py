#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_acute_HIA_pool_EPAregions.py
    python script to calculate asthma Hosp & asthma ED by region
    
    v2 - 01.04.21 - updated to look at acute outcomes by EPA region and season
    01.15.20 - added custom grid option
Created on Wed Sep 30 16:48:03 2020
@author: kodell
"""
#%% user inputs
# file locations
kPM_fp = '/Users/kodell/Desktop/data_not4gdrive/HIA_inputs/PM/'
kPM_desc = '_v2_medbk_final' # which kPM file to load
pop_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/population/rg_population_2010.nc'
br_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/acute_baseline_rates_US.csv'
state_grid_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/state_grid.npy'

# out file and figure path and descrption
out_fp = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/acute_outcomes/'
out_fig_path ='/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/acute_outcomes/'
out_desc = '_byBreyregion_kPM'+kPM_desc+'_final_reviewer_edits'

# years to calc HIA
years = [2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
       2017, 2018]

# set regions and assign states for custom option
region_names = ['Northwest','Midwest','Northeast',
                'Rocky Mountains','Great Plains','Mid Atlantic',
                'Southwest','Southern Plains','Southeast']
# for plotting with panel labels
region_names_panel = ['(a) Northwest, NW','(b) Midwest, MW','(c) Northeast, NE',
                '(d) Rocky Mountains, RM','(e) Great Plains, GP','(f) Mid Atlantic, MA',
                '(g) Southwest, SW','(h) Southern Plains, SP','(i) Southeast, SE']

# states for each region
region_states = [['WA','OR','ID'],
                 ['OH','IL','IN','MI','WI','MN'],
                 ['ME','NH','RI','VT','CT','MA','NJ','NY'],
                 ['MT','WY','UT','CO','NM'],
                 ['ND','SD','IA','KS','NE','MO'],
                 ['PA','WV','VA','MD','DE','DC'],
                 ['AZ','NV','CA'],
                 ['TX','AR','OK','LA'],
                 ['KY','TN','SC','NC','GA','AL','MS','FL']]

# assign region colors for plotting
region_colors = ['#1f78b4','#cab2d6','#a6cee3','#6a3d9a','#ffff99','#fdbf6f','#ff7f00','#fb9a99','#b2df8a']

#%% import modules
import numpy as np
import plotly.graph_objects as go
from netCDF4 import Dataset
import plotly.io as pio
from plotly.subplots import make_subplots
import pandas as pd
pio.renderers.default = "chrome"
pio.templates.default = "seaborn"

import matplotlib 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from wrf import (getvar, get_cartopy, cartopy_xlim,cartopy_ylim)

#%% user defined functions
# function to calculate HIAs
# odds ratio needs to be for delta PM of 10 ugm-3
# baseline rate needs to be daily
def calc_HIA(name,odds_ratio_lci,odds_ratio,odds_ratio_uci,baseline_rate,PM,pop):
    beta = np.log(odds_ratio)/10.0
    outcome = baseline_rate*(1.0-np.exp(-(beta)*PM))*pop
    # lci
    beta_lci = np.log(odds_ratio_lci)/10.0
    outcome_lci = baseline_rate*(1.0-np.exp(-(beta_lci)*PM))*pop
    # uci
    beta_uci = np.log(odds_ratio_uci)/10.0
    outcome_uci = baseline_rate*(1.0-np.exp(-(beta_uci)*PM))*pop

    return outcome_lci, outcome, outcome_uci

#%% load data
# PM data
kPM_fid = Dataset(kPM_fp + 'krigedPM25_06-18'+kPM_desc+'.nc')
allday_smokePM = kPM_fid['Smoke PM25'][:].data
lons = kPM_fid['lon'][:].data
lats = kPM_fid['lat'][:].data
kPM_fid.close()

# population data
nc_fid = Dataset(pop_file, 'r')
population = nc_fid.variables['population'][:].data
nc_fid.close()

# baseline rates
br = pd.read_csv(br_file)
br.set_index('outcome',inplace=True)

# state grid
state_grid = np.load(state_grid_file)

#%% make dataframe with EPI RRs
# use RRs from Borchers Arriagada et al., 2019 (10.1016/j.envres.2019.108777)
poolRRs = pd.DataFrame(data={'Study':['Borchers Arriagada et al., 2019','Borchers Arriagada et al., 2019'],
                             'Region':['USA meta-analysis','USA meta-analysis'],
                             'Outcome(s)':['asth ED','asth Hosp'],
                             'OR 10ug/m^3':[1.07,1.08],
                             'OR 95% LCI 10ug/m^3':[1.03,1.03],
                             'OR 95% UCI 10ug/m^3':[1.11,1.14]})

#%% set PM input for HIA
# test sensitivity to a 'threshold concentration' 
# allday_smokePM = np.where(allday_smokePM1 < 20, 0, allday_smokePM1)

pm_input = allday_smokePM # all US smoke, lag day zero

#%% set up region grid
num_regions = len(region_names)
epa_region_grid = np.copy(state_grid)
epa_region_grid[:] = 0
for ri in range(num_regions):
    states = region_states[ri]
    for stateID in states:
        inds=np.where(state_grid==stateID)
        epa_region_grid[inds] = ri+1

#%% mask values outside the US
non_US_inds = np.where(state_grid=='NA')
population[non_US_inds] = np.nan
allday_smokePM[:,:,non_US_inds[0],non_US_inds[1]] = np.nan

#%% calc HIA on all days and sum across the year
annual_event_all = np.empty([1,13,189,309])
annual_event_all_lci = np.empty([1,13,189,309])
annual_event_all_uci = np.empty([1,13,189,309])

JFM_all = np.empty([1,13,189,309])
AMJ_all = np.empty([1,13,189,309])
JAS_all = np.empty([1,13,189,309])
OND_all = np.empty([1,13,189,309])

outcome_names = []
for i in range(poolRRs.shape[0]):
    outcome_name = poolRRs.iloc[i]['Outcome(s)']
    
    # calculate via daily
    event_lci, event, event_uci = calc_HIA(poolRRs.iloc[i]['Study'],
                                           poolRRs.iloc[i]['OR 95% LCI 10ug/m^3'],
                                           poolRRs.iloc[i]['OR 10ug/m^3'],
                                           poolRRs.iloc[i]['OR 95% UCI 10ug/m^3'],
                                           br.loc[outcome_name]['baseline rate']/365.0,
                                           pm_input,
                                           population)
    outcome_names.append(outcome_name)
    
    # annual sums
    annual_event = np.array([np.nansum(event,axis=1)])
    annual_event_all = np.concatenate([annual_event_all,annual_event],axis=0)
   
    annual_event_lci = np.array([np.nansum(event_lci,axis=1)])
    annual_event_all_lci = np.concatenate([annual_event_all_lci,annual_event_lci],axis=0)
   
    annual_event_uci = np.array([np.nansum(event_uci,axis=1)])
    annual_event_all_uci = np.concatenate([annual_event_all_uci,annual_event_uci],axis=0)
    
    #seasonal sums
    # JFM, AMJ, JAS, OND
    JFM = np.array([np.nansum(event[:,:91,:,:],axis=1)])
    JFM_all = np.concatenate([JFM_all,JFM],axis=0)
    
    AMJ = np.array([np.nansum(event[:,91:182,:,:],axis=1)])
    AMJ_all = np.concatenate([AMJ_all,AMJ],axis=0)
    
    JAS = np.array([np.nansum(event[:,182:274,:,:],axis=1)])
    JAS_all = np.concatenate([JAS_all,JAS],axis=0)
    
    OND = np.array([np.nansum(event[:,274:,:,:],axis=1)])
    OND_all = np.concatenate([OND_all,OND],axis=0)
    '''
    # DJF, MAM, JJA, SON test in reponse to reviewers and new supplemental figures
    # going to use old abbr so we dont have to change that below
    JFM1 = np.array([np.nansum(event[:,:60,:,:],axis=1)])
    JFM2 = np.array([np.nansum(event[:,335:,:,:],axis=1)])
    JFM = JFM1 + JFM2
    JFM_all = np.concatenate([JFM_all,JFM],axis=0)
   
    AMJ = np.array([np.nansum(event[:,60:152,:,:],axis=1)])
    AMJ_all = np.concatenate([AMJ_all,AMJ],axis=0)
    
    JAS = np.array([np.nansum(event[:,152:244,:,:],axis=1)])
    JAS_all = np.concatenate([JAS_all,JAS],axis=0)
    
    OND = np.array([np.nansum(event[:,244:335,:,:],axis=1)])
    OND_all = np.concatenate([OND_all,OND],axis=0)
    '''
    i += 1

# sum across the grid cells, and remove first empty column
annual_event_all = annual_event_all[1:,:,:,:]
annual_event_all_lci = annual_event_all_lci[1:,:,:,:]
annual_event_all_uci = annual_event_all_uci[1:,:,:,:]

annual_event_sum = np.nansum(annual_event_all,axis=(2,3))
annual_event_sum_lci = np.nansum(annual_event_all_lci,axis=(2,3))
annual_event_sum_uci = np.nansum(annual_event_all_uci,axis=(2,3))

# stack seasons, and remove first empty column
byseason_all = np.stack([JFM_all[1:,:,:,:],AMJ_all[1:,:,:,:],
                         JAS_all[1:,:,:,:],OND_all[1:,:,:,:]],axis=2)

# for each outcome calculate total # of events
asth_hosp_total = np.nansum(population*br.loc['asth Hosp']['baseline rate'])
asth_ED_total = np.nansum(population*br.loc['asth ED']['baseline rate'])

#%% total events by EPA region
annual_event_sum_byEPAr = np.empty([annual_event_sum.shape[0],annual_event_sum.shape[1],num_regions])
annual_event_sum_byEPAr[:] = 0.0

byseason_sum_byEPAr = np.empty([annual_event_sum.shape[0],annual_event_sum.shape[1],4,num_regions])
byseason_sum_byEPAr[:] = 0.0

total_asth_hosp_byEPAr = np.empty([1,annual_event_sum.shape[1],num_regions])
total_asth_hosp_byEPAr[:] = 0.0

total_asth_ED_byEPAr = np.empty([1,annual_event_sum.shape[1],num_regions])
total_asth_ED_byEPAr[:] = 0.0

for i in range(num_regions):
    region_inds = np.where(epa_region_grid==str(i+1))
    annual_event_sum_byEPAr[:,:,i] = np.nansum(annual_event_all[:,:,region_inds[0],region_inds[1]],axis=(2))
    byseason_sum_byEPAr[:,:,:,i] = np.nansum(byseason_all[:,:,:,region_inds[0],region_inds[1]],axis=(3))
    region_pop = np.nansum(population[region_inds[0],region_inds[1]])
    total_asth_hosp_byEPAr[:,:,i] = region_pop*br.loc['asth Hosp']['baseline rate']
    total_asth_ED_byEPAr[:,:,i] = region_pop*br.loc['asth ED']['baseline rate']

#%% print numbers that are in paper
# total asthma ED visits, hospital admissions
print('tot ED by year, range',annual_event_sum[0,:].min(),annual_event_sum[0,:].max())
print('tot hosp by year, range',annual_event_sum[1,:].min(),annual_event_sum[1,:].max())
# percent att. to smoke
print('% ED by year, range',100.0*(annual_event_sum[0,:].min()/asth_ED_total),
      100.0*(annual_event_sum[0,:].max()/asth_ED_total))
print('% hosp by year,range',100.0*(annual_event_sum[1,:].min()/asth_hosp_total),
      100.0*(annual_event_sum[1,:].max()/asth_hosp_total))
'''
print('tot ED by year, lci range',annual_event_sum_lci[0,:].min(),annual_event_sum_lci[0,:].max())
print('tot hosp by year, lci range',annual_event_sum_lci[1,:].min(),annual_event_sum_lci[1,:].max())

print('tot ED by year, uci range',annual_event_sum_uci[0,:].min(),annual_event_sum_uci[0,:].max())
print('tot hosp by year, uci range',annual_event_sum_uci[1,:].min(),annual_event_sum_uci[1,:].max())
'''

# eUS asthma ED visits, hospital admissions
region_names=np.array(region_names)
eUS_rinds = []
for rname in ['Southeast','Midwest','Northeast','Great Plains',
                                      'Mid Atlantic','Southern Plains']:
    eUS_rinds.append(np.where(region_names==rname)[0][0])
    
eUS_annual_event_sum = np.sum(annual_event_sum_byEPAr[:,:,eUS_rinds],axis=2)
eUS_asth_ED_total = np.sum(total_asth_ED_byEPAr[:,:,eUS_rinds],axis=2)[0][0]
eUS_asth_hosp_total = np.sum(total_asth_hosp_byEPAr[:,:,eUS_rinds],axis=2)[0][0]
print('eUS tot ED by year, avg',eUS_annual_event_sum[0,:].mean())
print('eUS tot hosp by year, avg',eUS_annual_event_sum[1,:].mean())
# eUS percent att. to smoke
print('eUS % ED by year, avg',100.0*(eUS_annual_event_sum[0,:].mean()/eUS_asth_ED_total))
print('eUS % hosp by year, avg',100.0*(eUS_annual_event_sum[1,:].mean()/eUS_asth_hosp_total))

# wUS asthma ED visits, hospital admissions
wUS_rinds = []
for rname in ['Southwest','Rocky Mountains','Northwest']:
    wUS_rinds.append(np.where(region_names==rname)[0][0])
    
wUS_annual_event_sum = np.sum(annual_event_sum_byEPAr[:,:,wUS_rinds],axis=2)
wUS_asth_ED_total = np.sum(total_asth_ED_byEPAr[:,:,wUS_rinds],axis=2)[0][0]
wUS_asth_hosp_total = np.sum(total_asth_hosp_byEPAr[:,:,wUS_rinds],axis=2)[0][0]
print('wUS tot ED by year, avg',wUS_annual_event_sum[0,:].mean())
print('wUS tot hosp by year, avg',wUS_annual_event_sum[1,:].mean())
# wUS percent att. to smoke
print('wUS % ED by year, avg',100.0*(wUS_annual_event_sum[0,:].mean()/wUS_asth_ED_total))
print('wUS % hosp by year, avg',100.0*(wUS_annual_event_sum[1,:].mean()/wUS_asth_hosp_total))

#%% save data
out_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/acute_outcomes/acute_events_'+out_desc
np.savez(out_fn,annual_event_all=annual_event_all,byseason_all=byseason_all,
         lons=lons,lats=lats)
print('data saved, making figures')

#%% make figures
# 1) map of regions, included in Figure 2 
# (region abbrivations added in powerpoint)
cmap = matplotlib.colors.ListedColormap(region_colors)
cmap.set_under('white')
bounds = np.arange(1,num_regions+2)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

ncfile = Dataset('/Users/kodell/Local Google Drive /wrfout_d01_2016-07-22_00:00:00')
slp = getvar(ncfile,"slp")
cart_proj = get_cartopy(slp)
ncfile.close()

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection=cart_proj)
ax.set_extent([-125, -66.5, 20, 50], ccrs.Geodetic())
plt.gca().outline_patch.set_visible(False)
ax.patch.set_visible(False)
shapename = 'admin_1_states_provinces_lakes_shp'
states_shp = shpreader.natural_earth(resolution='110m',category='cultural', name=shapename)
facecolor = (0, 0, 0, 0)
edgecolor = 'black'
for state in shpreader.Reader(states_shp).records():
    sID = state.attributes['postal']
    if sID in ['AK','HI']:
        continue
    for i in range(len(region_states)):
        if sID in region_states[i]:
            facecolor = region_colors[i]
    ax.add_geometries([state.geometry], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)
plt.savefig(out_fig_path+'region_map'+out_desc+'_paper.png',dpi=300)
plt.show()

#%% 2) events by EPA region annually
# Figure 2 and S3
outcome_names_fig = ['Asthma ED vists', 'Asthma hospital<br>admissions'] 
outcome_names_save = ['Asthma_ED_Vists', 'Asthma_Hospitalizations']
fig = make_subplots(rows=2,cols=1,specs=[[{"secondary_y": True}],[{"secondary_y": True}]])

# remove negatives before plotting
annual_event_sum_byEPAr_plot = np.where(annual_event_sum_byEPAr<0,0,annual_event_sum_byEPAr)
for j in [1,2,4,5,7,8,0,3,6]: # to organize order on plot
        fig.add_trace(go.Bar(x=years,y=annual_event_sum_byEPAr_plot[0,:,j],
                    name=region_names[j],
                     marker_color = region_colors[j],showlegend=False),
                     row=1,col=1,secondary_y=False)
        fig.add_trace(go.Bar(x=years,y=annual_event_sum_byEPAr_plot[1,:,j],
                    name=region_names[j],
                     marker_color = region_colors[j],showlegend=False),
                     row=2,col=1,secondary_y=False)
        fig.add_trace(go.Bar(x=years,y=100.0*(annual_event_sum_byEPAr_plot[0,:,j]/asth_ED_total),
                    name=region_names[j],
                     marker_color = region_colors[j],showlegend=False),
                     row=1,col=1,secondary_y=True)
        fig.add_trace(go.Bar(x=years,y=100.0*(annual_event_sum_byEPAr_plot[1,:,j]/asth_hosp_total),
                    name=region_names[j],
                     marker_color = region_colors[j],showlegend=False),
                     row=2,col=1,secondary_y=True)
# add error bars to both y-axes for bar allignment
fig.add_trace(go.Bar(x=years,y=[0]*len(years),marker_color = 'white',showlegend=False,
                     error_y=dict(type='data',color='dimgrey',
                                  array=np.sum(annual_event_sum_byEPAr_plot[0,:,:],axis=1)-annual_event_sum_lci[0,:])),
              row=1,col=1,secondary_y=False)
fig.add_trace(go.Bar(x=years,y=[0]*len(years),marker_color = 'white',showlegend=False,
                     error_y=dict(type='data',color='dimgrey',
                                  array=np.sum(annual_event_sum_byEPAr_plot[1,:,:],axis=1)-annual_event_sum_lci[1,:])),
              row=2,col=1,secondary_y=False)
fig.add_trace(go.Bar(x=years,y=[0]*len(years),marker_color = 'white',showlegend=False,
                     error_y=dict(type='data',color='dimgrey',
                                  array=100.0*(np.sum(annual_event_sum_byEPAr_plot[0,:,:],axis=1)-annual_event_sum_lci[0,:])/asth_ED_total)),
                     row=1,col=1,secondary_y=True)
fig.add_trace(go.Bar(x=years,y=[0]*len(years),marker_color = 'white',showlegend=False,
                     error_y=dict(type='data',color='dimgrey',
                                  array=100.0*(np.sum(annual_event_sum_byEPAr_plot[1,:,:],axis=1)-annual_event_sum_lci[1,:])/asth_hosp_total)),
                     row=2,col=1,secondary_y=True)
fig.update_layout(plot_bgcolor='white',barmode='stack',font_family='Arial',font_size=18,
                  title='<b>2006-2018 smoke-attributable asthma morbidity by region</b>')
fig.update_yaxes(title_text=outcome_names_fig[0],ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 showline=True, linewidth=1, linecolor='black',
                 row=1,col=1,secondary_y=False)
fig.update_yaxes(title_text=outcome_names_fig[1],ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 showline=True, linewidth=1, linecolor='black',
                 row=2,col=1,secondary_y=False)
fig.update_yaxes(title_text='% of total asthma<br>ED visits',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 showline=True, linewidth=1, linecolor='black',
                 row=1,col=1,secondary_y=True)
fig.update_yaxes(title_text='% of total asthma<br> hospital admissions',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 showline=True, linewidth=1, linecolor='black',
                 row=2,col=1,secondary_y=True)
fig.update_xaxes(showline=True, linewidth=1, linecolor='black',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10)

fig.show()
fig.write_image(out_fig_path + 'annual_acute_'+out_desc+'.png',scale=4,height=600,width=1000)#width=800)

#%% 3a) Events by EPA region by season
# not in paper
# before plotting, remove negatives ... plotly is weird about these
plot_byseason_sum_byEPAr = np.where(byseason_sum_byEPAr<0,0,byseason_sum_byEPAr)
fig = make_subplots(rows=2,cols=2,subplot_titles = ['Winter (JFM)','Spring (AMJ)','Summer (JAS)','Fall (OND)'])
for j in [1,2,4,5,7,8,0,3,6]:
        fig.add_trace(go.Bar(x=years,y=plot_byseason_sum_byEPAr[0,:,0,j],name=region_names[j],
                     marker_color = region_colors[j],showlegend=False),
                     row=1,col=1)
        fig.add_trace(go.Bar(x=years,y=plot_byseason_sum_byEPAr[0,:,1,j],name=region_names[j],
                     marker_color = region_colors[j],showlegend=False),
                     row=1,col=2)
        fig.add_trace(go.Bar(x=years,y=plot_byseason_sum_byEPAr[0,:,2,j],name=region_names[j],
                     marker_color = region_colors[j],showlegend=False),
                     row=2,col=1)
        fig.add_trace(go.Bar(x=years,y=plot_byseason_sum_byEPAr[0,:,3,j],name=region_names[j],
                     marker_color = region_colors[j],showlegend=False),
                     row=2,col=2)

fig.update_layout(plot_bgcolor='white',barmode='stack',
                  title = '<b>Smoke-attributable asthma ED vistis by season</b>',font_family='Arial')
fig.update_yaxes(title=outcome_names_fig[0],ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 showline=True, linewidth=1, linecolor='black')
fig.update_xaxes(ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 showline=True, linewidth=1, linecolor='black')
fig.write_image(out_fig_path + outcome_names_save[0]+out_desc+'_byseason.png',scale=4,height=600,width=800)
fig.show()

# 3b) plotted as a percentage of all events in that region
# not in paper
inds = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
show_leg = [True,False,False,False,False,False,False,False,False]
fig = make_subplots(rows=3,cols=3,subplot_titles = region_names_panel)
i = 0
for j in range(num_regions):
        fig.add_trace(go.Bar(x=years,y=100.0*(plot_byseason_sum_byEPAr[i,:,0,j]/np.nansum(plot_byseason_sum_byEPAr[i,:,:,j],axis=1)),
                     name='Winter (JFM)',
                     marker_color = 'blue',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1])
        fig.add_trace(go.Bar(x=years,y=100.0*(plot_byseason_sum_byEPAr[i,:,1,j]/np.nansum(plot_byseason_sum_byEPAr[i,:,:,j],axis=1)),
                     name='Spring (AMJ)',
                     marker_color = 'pink',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1])
        fig.add_trace(go.Bar(x=years,y=100.0*(plot_byseason_sum_byEPAr[i,:,2,j]/np.nansum(plot_byseason_sum_byEPAr[i,:,:,j],axis=1)),
                     name='Summer (JAS)',
                     marker_color = 'green',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1])
        fig.add_trace(go.Bar(x=years,y=100.0*(plot_byseason_sum_byEPAr[i,:,3,j]/np.nansum(plot_byseason_sum_byEPAr[i,:,:,j],axis=1)), 
                     name='Fall (OND)',
                     marker_color = 'orange',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1])
fig.update_layout(plot_bgcolor='white',barmode='stack',font_family='Arial',legend_traceorder='normal',
                  title = '<b>Normalized regional smoke-attributable asthma morbidity by season</b>')
fig.update_yaxes(showline=True, linewidth=1, linecolor='black',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10)
fig.update_yaxes(title='% of smoke-attributable<br>asthma ED visits',row=1,col=1)
fig.update_yaxes(title='% of smoke-attributable<br>asthma ED visits',row=2,col=1)
fig.update_yaxes(title='% of smoke-attributable<br>asthma ED visits',row=3,col=1)
fig.update_xaxes(showline=True, linewidth=1, linecolor='black',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10)
fig.show()
fig.write_image(out_fig_path + outcome_names_save[i]+out_desc+'_byseason2.png',scale=4,height=600,width=800)

#%% 3c) regional events normalized by events per year per region, colored by season
# Figure 3, and S4
inds = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
show_leg = [False,True,False,False,False,False,False,False,False]
#show_leg = [True,True,True,True,True,True,True,True,True,True]
fig = make_subplots(rows=3,cols=3,subplot_titles = region_names_panel,
                    specs=[[{"secondary_y": True},{"secondary_y": True},{"secondary_y": True}],
                           [{"secondary_y": True},{"secondary_y": True},{"secondary_y": True}],
                           [{"secondary_y": True},{"secondary_y": True},{"secondary_y": True}]])
for j in range(num_regions):
        # first plot number
        fig.add_trace(go.Bar(x=years,y=plot_byseason_sum_byEPAr[i,:,0,j],
                             name='Winter (JFM)',
                     marker_color = 'blue',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1],secondary_y=False)
        fig.add_trace(go.Bar(x=years,y=plot_byseason_sum_byEPAr[i,:,1,j],
                             name='Sping (AMJ)',
                     marker_color = 'pink',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1],secondary_y=False)
        fig.add_trace(go.Bar(x=years,y=plot_byseason_sum_byEPAr[i,:,2,j],
                             name='Summer (JAS)',
                     marker_color = 'green',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1],secondary_y=False)
        fig.add_trace(go.Bar(x=years,y=plot_byseason_sum_byEPAr[i,:,3,j],
                             name='Fall (OND)',
                     marker_color = 'orange',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1],secondary_y=False)
        
        # then overlay % on secondary y (always have to check with and without overlay)
        fig.add_trace(go.Bar(x=years,y=100.0*(plot_byseason_sum_byEPAr[i,:,0,j]/total_asth_ED_byEPAr[0,:,j]),name='Winter (JFM)',
                     marker_color = 'blue',showlegend=False),
                     row=inds[j][0],col=inds[j][1],secondary_y=True)
        fig.add_trace(go.Bar(x=years,y=100.0*(plot_byseason_sum_byEPAr[i,:,1,j]/total_asth_ED_byEPAr[0,:,j]),name='Sping (AMJ)',
                     marker_color = 'pink',showlegend=False),
                     row=inds[j][0],col=inds[j][1],secondary_y=True)
        fig.add_trace(go.Bar(x=years,y=100.0*(plot_byseason_sum_byEPAr[i,:,2,j]/total_asth_ED_byEPAr[0,:,j]),name='Summer (JAS)',
                     marker_color = 'green',showlegend=False),
                     row=inds[j][0],col=inds[j][1],secondary_y=True)
        fig.add_trace(go.Bar(x=years,y=100.0*(plot_byseason_sum_byEPAr[i,:,3,j]/total_asth_ED_byEPAr[0,:,j]),name='Fall (OND)',
                     marker_color = 'orange',showlegend=False),
                     row=inds[j][0],col=inds[j][1],secondary_y=True)
fig.update_layout(plot_bgcolor='white',barmode='stack',font_family='Arial',legend_traceorder='normal',
                  title='<b>Regional smoke-attributable asthma morbidity by season</b>')
fig.update_yaxes(showline=True, linewidth=1, linecolor='black',ticks="outside", tickwidth=2, tickcolor='gray', 
                 ticklen=10, secondary_y=False)
fig.update_yaxes(showline=True, linewidth=1, linecolor='crimson',ticks="outside", tickwidth=2, tickcolor='crimson', 
                 ticklen=10, tickfont=dict(family='Arial', color='crimson'), secondary_y=True)
fig.update_yaxes(title=outcome_names_fig[i],row=1,col=1,secondary_y=False)
fig.update_yaxes(title=outcome_names_fig[i],row=2,col=1,secondary_y=False)
fig.update_yaxes(title=outcome_names_fig[i],row=3,col=1,secondary_y=False)

fig.update_yaxes(title='% of region total<br>asthma ED visits',title_font=dict(color='crimson',family = 'Arial'),
                 linecolor='red',row=1,col=3,secondary_y=True)
fig.update_yaxes(title='% of region total<br>asthma ED visits',title_font=dict(color='crimson',family='Arial'),
                 linecolor='red',row=2,col=3,secondary_y=True)
fig.update_yaxes(title='% of region total<br>asthma ED visits',title_font=dict(color='crimson',family='Arial'),
                 linecolor='red',row=3,col=3,secondary_y=True)

fig.update_xaxes(showline=True, linewidth=1, linecolor='black',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10,
                 tickvals=[2006,2010,2014,2018])
fig.update_layout(legend=dict(orientation="v",yanchor="top",y=1,xanchor="left",x=1.1))
fig.show()
fig.write_image(out_fig_path + outcome_names_save[i]+out_desc+'_byseason2.png',scale=4,height=600,width=900)

#%% 3d) plotted as a percentage of all events in that region for entire time period
# Figure S5
fig = make_subplots(rows=1,cols=1)
sum_allyears = np.nansum(plot_byseason_sum_byEPAr[0,:,:,:],axis=0)
fig.add_trace(go.Bar(x=region_names,y=100.0*(sum_allyears[0,:]/np.nansum(sum_allyears,axis=0)),
             name='Winter (JFM)',
             marker_color = 'blue',showlegend=show_leg[j]))
fig.add_trace(go.Bar(x=region_names,y=100.0*(sum_allyears[1,:]/np.nansum(sum_allyears,axis=0)),
             name='Spring (AMJ)',
             marker_color = 'pink',showlegend=show_leg[j]))
fig.add_trace(go.Bar(x=region_names,y=100.0*(sum_allyears[2,:]/np.nansum(sum_allyears,axis=0)),
             name='Summer (JAS)',
             marker_color = 'green',showlegend=show_leg[j]))
fig.add_trace(go.Bar(x=region_names,y=100.0*(sum_allyears[3,:]/np.nansum(sum_allyears,axis=0)), 
             name='Fall (OND)',
             marker_color = 'orange',showlegend=show_leg[j]))
fig.update_layout(plot_bgcolor='white',barmode='stack',font_family='Arial',legend_traceorder='normal',
                  title = '<b>Normalized regional smoke-attributable asthma morbidity by season 2006-2018</b>')
fig.update_yaxes(showline=True, linewidth=1, linecolor='black',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10)
fig.update_yaxes(title='% of smoke-attributable<br>asthma ED visits')
fig.update_xaxes(showline=True, linewidth=1, linecolor='black',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10)
fig.show()
fig.write_image(out_fig_path + outcome_names_save[i]+out_desc+'_byseason_allyears.png',scale=4,height=600,width=800)

#%% plot graveyard
'''

# 4) time series of EDs for each year - just so we know what this looks like 
# a little noisy - but sure enough dual peak in spring and summer in general
fig = go.Figure()
for i in range(len(years)):
    fig.add_trace(go.Scatter(x=np.arange(0,366),y=np.nansum(event[i,:,:,:],axis=(1,2)),
                             mode = 'lines',name=str(years[i])))
fig.update_layout(title='timeseries of asthma hosp admissions')
fig.show()


fig.add_trace(go.Bar(x=years,y=wUS_annual_event_sum[i,:],name=outcome_names_fig[i],
                     marker_color = colors[i],
                     width = 0.3,offset = offsets[i],base=0))
fig.add_trace(go.Bar(x=years,y=eUS_annual_event_sum[i,:],name=outcome_names_fig[i],
                     marker_color = colors[i],showlegend=False,opacity=0.7,
                     width = 0.3,offset = offsets[i],base=wUS_annual_event_sum[i,:],
                     error_y=dict(type='data',color='dimgrey',
                                  array=annual_event_sum[i,:]-annual_event_sum_lci[i,:])))


fig = go.Figure()
colors = ['indianred','cornflowerblue']#,'mediumaquamarine']
#for i in range(3):
for i in range(2):
    fig.add_trace(go.Bar(x=years,y=annual_event_sum[i,:],
                         name=outcome_names_fig[i],marker_color = colors[i],
                         error_y=dict(type='data',color='dimgrey',
                                      array=annual_event_sum[i,:]-annual_event_sum_lci[i,:])))
fig.update_layout(plot_bgcolor='white')
fig.update_yaxes(title='annual events')
fig.show()
fig.write_image(out_fig_path + 'acute_HIA_totals'+out_desc+'_total.png',scale=4,height=600,width=800)


fig = go.Figure()
totals = [asth_ED_total,asth_hosp_total]#,COPD_ED_total]
#for i in range(3):
for i in range(2):
    fig.add_trace(go.Bar(x=years,y=100.0*(annual_event_sum[i,:]/totals[i]),
                         name=outcome_names_fig[i],marker_color = colors[i],
                         error_y=dict(type='data',color='dimgrey',
                                      array=100.0*((annual_event_sum[i,:]/totals[i])-(annual_event_sum_lci[i,:]/totals[i])))))
fig.update_layout(plot_bgcolor='white')
fig.update_yaxes(title='percent of annual events<br>attributable to smoke PM<sub>2.5<sub>')
fig.show()
fig.write_image(out_fig_path + 'acute_HIA_fractions'+out_desc+'.png',scale=4,height=600,width=800)

# 3) fraction of total events by season
inds = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
show_leg = [True,False,False,False,False,False,False,False,False]
fig = make_subplots(rows=3,cols=3,subplot_titles = region_names)
i = 0
for j in range(num_regions):
        fig.add_trace(go.Bar(x=years,y=byseason_sum_byEPAr[i,:,0,j],name='JFM',
                     marker_color = 'blue',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1])
        fig.add_trace(go.Bar(x=years,y=byseason_sum_byEPAr[i,:,1,j],name='AMJ',
                     marker_color = 'pink',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1])
        fig.add_trace(go.Bar(x=years,y=byseason_sum_byEPAr[i,:,2,j],name='JAS',
                     marker_color = 'green',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1])
        fig.add_trace(go.Bar(x=years,y=byseason_sum_byEPAr[i,:,3,j],name='OND',
                     marker_color = 'orange',showlegend=show_leg[j]),
                     row=inds[j][0],col=inds[j][1])
fig.update_layout(plot_bgcolor='white',barmode='stack',font_family='Tahoma',legend_traceorder='normal')
fig.update_yaxes(showline=True, linewidth=1, linecolor='black',ticks="outside", tickwidth=2, tickcolor='gray', ticklen=10)
fig.update_yaxes(title=outcome_names_fig[i],row=1,col=1)
fig.update_yaxes(title=outcome_names_fig[i],row=2,col=1)
fig.update_yaxes(title=outcome_names_fig[i],row=3,col=1)
fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
fig.show()
fig.write_image(out_fig_path + outcome_names_save[i]+out_desc+'_byseason2.png',scale=4,height=600,width=800)

# plot to check regions
cmap = matplotlib.colors.ListedColormap(region_colors)
cmap.set_under('white')
bounds = np.arange(1,num_regions+2)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

ncfile = Dataset('/Users/kodell/Local Google Drive /wrfout_d01_2016-07-22_00:00:00')
slp = getvar(ncfile,"slp")
cart_proj = get_cartopy(slp)

fig = plt.figure()
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
cs = ax.pcolor(lons,lats, epa_region_grid,cmap=cmap,norm=norm,shading='nearest',
               transform = ccrs.PlateCarree())

plt.title('regions',fontsize=14)
plt.savefig(out_fig_path+'region_map'+out_desc+'.png',dpi=300)
plt.show()

'''

    