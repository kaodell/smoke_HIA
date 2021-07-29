# Calculate mortalities
# Bilsback 2019
# redone to plot with plotly, by Kate O'Dell, June 2020
# updated to calc mortality by state as well (originally in calc_by_country)
#%% import modules
import numpy as np
import netCDF4
import plotly
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "chrome"
pio.templates.default = "seaborn"
import cartopy.feature as cfeature
import pandas as pd
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.feature import NaturalEarthFeature
from wrf import (getvar, get_cartopy, cartopy_xlim,cartopy_ylim)
import matplotlib.pyplot as plt
import matplotlib

# set figure fonts
from matplotlib.font_manager import FontProperties

font = FontProperties()
font.set_family('sansserif')
font.set_name('Tahoma')
#%% user inputs
desc = '_v2_medbk_final_2010pop'
pop_file = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/population/rg_population_2010.nc'
id_grid_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/inputs/state_grid.npy'
mort_by_grid_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/mortalities/wilfire_smoke_premature_mort_gemm_'+desc+'.nc'

# load all cause or 5 leading causes?
cause = 'allcause'
# cause = '5leadingcauses'

# description for output files
out_desc = desc + cause + '_reviewer_edits'
# out figure path
out_data_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/mortalities/'
out_fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/mort_bystate/'

#%% user-defined functions
def plot_background(ax):
    ax.set_extent([235., 290., 20., 55.])
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    return ax
#%% load data
# load id grid
id_grid = np.load(id_grid_fn)

# load by-grid calculated data
nc_fid = netCDF4.Dataset(mort_by_grid_fn)
lon = nc_fid.variables['lon'][:]
lat = nc_fid.variables['lat'][:]
if cause == 'allcause':
    allPM = nc_fid.variables['total_mort_ac'][1,:,:]
    nosmokePM = nc_fid.variables['nosmoke_mort_ac'][1,:,:]
    smokePMaf = nc_fid.variables['smoke_mort_att_fraction_ac'][1,:,:]
if cause == '5leadingcauses':
    allPM = nc_fid.variables['total_mort'][1,:,:]
    nosmokePM = nc_fid.variables['nosmoke_mort'][1,:,:]
    smokePMaf = nc_fid.variables['smoke_mort_att_fraction'][1,:,:]
nc_fid.close()

# load population so we can sum by state
nc_pop = netCDF4.Dataset(pop_file, 'r')
population = nc_pop.variables['population'][:]

#%% sum deaths by state 
state_codes = []
allPM_state = []
nosmokePM_state = []
smokePMaf_state = []
state_pop = []
CDC_pops = []

gstate_allPM_mort = np.empty([189,309])
gstate_allPM_mort[:] = np.nan
gstate_smokePMaf_mort = np.empty([189,309])
gstate_smokePMaf_mort[:] = np.nan
gstate_pop = np.empty([189,309])
gstate_pop[:] = np.nan
wmort = 0
emort = 0
wpop = 0
epop = 0
for state_ID in np.unique(id_grid):
    print(state_ID)
    if state_ID=='NA':
        continue
    # calculate total by state
    allPM_total = np.nansum(np.where(id_grid == state_ID, allPM, 0))
    nosmokePM_total = np.nansum(np.where(id_grid == state_ID, nosmokePM, 0))
    smokePMaf_total = np.nansum(np.where(id_grid == state_ID, smokePMaf, 0))
    pop_total = np.nansum(np.where(id_grid == state_ID, population, 0))

    # add to arrays
    allPM_state.append(allPM_total)
    nosmokePM_state.append(nosmokePM_total)
    smokePMaf_state.append(smokePMaf_total)
    state_pop.append(pop_total)
    state_codes.append(state_ID)    
    
    if state_ID in ['CA','WA','OR','MT','ID','WY','CO','UT','NV','AZ','NM']:
        wmort += smokePMaf_total
        wpop += pop_total
    else:
        emort += smokePMaf_total
        epop += pop_total

    # add to grid for plotting for paper
    state_inds = np.where(id_grid==state_ID)
    gstate_allPM_mort[state_inds] = allPM_total
    gstate_smokePMaf_mort[state_inds] = smokePMaf_total
    gstate_pop[state_inds] = pop_total
    
# make sure all deaths in contig US are being added to a state
print('total US smoke deaths should match: ',np.sum(smokePMaf_state),np.nansum(smokePMaf))
print('total east mort',emort,'total west',wmort)
print('pct east mort',emort/(epop*0.0080095),'pct west mort',wmort/(wpop*0.0080095))
#%% calc percent deaths
# by all-cause mortality rate 800.95/100,000
# let's do deaths per 100,000 in the state
allPM_state_popw = 100000*np.array(allPM_state)/(np.array(state_pop))
nosmokePM_state_popw = 100000*np.array(nosmokePM_state)/(np.array(state_pop))
smokePMaf_state_popw = 100000*np.array(smokePMaf_state)/(np.array(state_pop))

# all-cause death rate is 800.95 per 100,000 from http://ghdx.healthdata.org/gbd-results-tool (2019 analysis)
# on webpage select: location: Location: US, Year: 2010 Metric: Rate, Age: All ages, Measure: Death, Sex: Both
# note this is different than 'all cause' above which is non-communicable diseases + LRI (defined in Burnett et al., 2018)
allPM_state_mortw = 100.0*np.array(allPM_state)/(np.array(state_pop)*0.0080095)
nosmokePM_state_mortw = 100.0*np.array(nosmokePM_state)/(np.array(state_pop)*0.0080095)
smokePMaf_state_mortw = 100.0*np.array(smokePMaf_state)/(np.array(state_pop)*0.0080095)

# now for gridded version (this is what we plot for the paper)
gstate_allPM_mortw = 100.0*gstate_allPM_mort/(gstate_pop*0.0080095)
gstate_smokePMaf_mortw = 100.0*gstate_smokePMaf_mort/(gstate_pop*0.0080095)

#%% save data
state_totals_df = pd.DataFrame(data = {'state':state_codes,
                                       'population':state_pop,
                                       'total PM deaths':allPM_state,
                                       'no smoke PM deaths':nosmokePM_state,
                                       'smoke PM deaths (att. fraction)':smokePMaf_state,
                                       'percent deaths smoke':smokePMaf_state_mortw})
state_totals_df.to_csv(out_data_path + 'mort_by_state' + out_desc + '.csv')

#%% Plot mortalities plotly
fig = go.Figure(data=go.Choropleth(
    locations=state_codes, # Spatial coordinates
    z = allPM_state, # Data to be color-coded
    locationmode = 'USA-states', # set of locations match entries in `locations`
    colorscale = 'Reds',
    colorbar_title = "mortalities per year"))

# all PM mortalities
fig.update_layout(title_text = 'all PM mortalities',geo_scope='usa',) # limit map scope to USA)
fig.show()
fig.write_image(out_fig_path+'allPMmort_bystate_'+out_desc+'.png',scale=8)
# mortalities not attributable to smoke
fig = go.Figure(data=go.Choropleth(locations=state_codes, z = nosmokePM_state, 
    locationmode = 'USA-states', colorscale = 'Reds',colorbar_title = "mortalities per year"))
fig.update_layout(title_text = 'nosmoke PM mortalities',geo_scope='usa',)
fig.show()
fig.write_image(out_fig_path+'nosmokePMmort_bystate_'+out_desc+'.png',scale=8)  
# smoke attributable fraction
fig = go.Figure(data=go.Choropleth(locations=state_codes, z = smokePMaf_state, 
    locationmode = 'USA-states', colorscale = 'Reds',colorbar_title = "mortalities per year"))
fig.update_layout(title_text = 'smoke PM mortalities',geo_scope='usa',) 
fig.show()
fig.write_image(out_fig_path+'smokePMmort_bystate_'+out_desc+'.png',scale=8)  

# population weighted versions of figure
fig = go.Figure(data=go.Choropleth(locations=state_codes, z = allPM_state_mortw,
    locationmode = 'USA-states',colorscale = 'Reds',
    colorbar_title = "% of mortalities per year",zmin=4,zmax=12))
fig.update_layout(title_text = 'all PM mortalities',geo_scope='usa',) 
fig.show()
fig.write_image(out_fig_path+'allPMmort_mortweight_bystate_'+out_desc+'.png',scale=8)  

fig = go.Figure(data=go.Choropleth(locations=state_codes,z = nosmokePM_state_mortw , 
    locationmode = 'USA-states', colorscale = 'Reds',
    colorbar_title = "% of mortalities per year",zmin=4,zmax=12))
fig.update_layout(title_text = 'nosmoke PM mortalities',geo_scope='usa',)
fig.show()
fig.write_image(out_fig_path+'nosmokePMmort_mortweight_bystate_'+out_desc+'.png',scale=8)  

fig = go.Figure(data=go.Choropleth(locations=state_codes, z = smokePMaf_state_mortw,
    locationmode = 'USA-states', colorscale = 'Reds',
    colorbar_title = "% of mortalities per year",zmin=0,zmax=1.5))
fig.update_layout(title_text = 'smoke PM mortalities',geo_scope='usa',) 
fig.show()
fig.write_image(out_fig_path+'smokePMmort_mortweight_bystate_'+out_desc+'.png',
                scale=8)  


#%% make plot for the paper

# wrf plotting info: https://wrf-python.readthedocs.io/en/latest/plot.html 
# load wrf file to get cartopy projection
ncfile = netCDF4.Dataset('/Users/kodell/Local Google Drive /wrfout_d01_2016-07-22_00:00:00')
slp = getvar(ncfile,"slp")
cart_proj = get_cartopy(slp)

def mk_map(ax):
    ax.patch.set_visible(False)
    shapename = 'admin_1_states_provinces_lakes_shp'
    states_shp = shpreader.natural_earth(resolution='110m',category='cultural', name=shapename)
    facecolor = (0, 0, 0, 0)
    edgecolor = 'black'
    for state in shpreader.Reader(states_shp).geometries():
        ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)


# set up a blank map with multiple subplots
fig,axarr = plt.subplots(ncols=2,nrows=2,figsize=(5,5),
                      subplot_kw={'projection': cart_proj})
titles = np.array([['(a) Total PM$_{2.5}$\n mortalities','(b) Smoke PM$_{2.5}$\n attributable mortalities'],
                   ['(c) Percent of all mortalities\n attributable to total PM$_{2.5}$',
                    '(d) Percent of all mortalities\n attributable to smoke PM$_{2.5}$']])
data = np.array([[gstate_allPM_mort/1000.0,gstate_smokePMaf_mort/1000.0],
                   [gstate_allPM_mortw,gstate_smokePMaf_mortw]])
cbar_labels = np.array([['Annual mortalities (1000s)','Annual mortalities (1000s)'],
                   ['Percent annual mortalities (%)','Percent annual mortalities (%)']])
cmin = np.array([[1,0.01],
                   [6.0,0.1]])
cmax = np.array([[33,0.85],
                   [11.0,1.2]])
cmaps = np.array([['Purples','Oranges'],['Purples','Oranges']])
for i in range(2):
    for j in range(2):
        ax = axarr[i,j]
        mk_map(ax)
        ax.outline_patch.set_edgecolor('white')
        cs = ax.pcolor(lon,lat, data[i,j], transform=ccrs.PlateCarree(),shading='nearest',cmap='Reds',
                       vmin=cmin[i,j],vmax=cmax[i,j])
        cax,kw = matplotlib.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.7)
        cbar=fig.colorbar(cs,cax=cax,extend='min',**kw)
        cbar.set_label(cbar_labels[i,j],fontsize=8,fontproperties = font)
        ax.set_title(titles[i,j],fontsize=10,fontproperties = font)
plt.subplots_adjust(top=1.0,
bottom=0.17,
left=0.12,
right=0.89,
hspace=0.0,
wspace=0.2)
plt.savefig(out_fig_path+'paper_mort_map'+out_desc+'_updateshp.png',dpi=300)
plt.show()

 
#%% plots to check clororpleth
# use catopy to plot arrays where total is assigned to each grid cell
# looks great!
'''
fig, axarr = plt.subplots(nrows=3, ncols=1, figsize=(20, 13), constrained_layout=True,
                      subplot_kw={'projection': ccrs.PlateCarree()})
axlist = axarr.flatten()
for ax in axlist:
    plot_background(ax)
c=axlist[0].contourf(lon, lat, allPM,60, transform=ccrs.PlateCarree())
axlist[0].set_title('total PM mort', fontsize=16)
cb1 = fig.colorbar(c, ax=axlist[0], orientation='horizontal', shrink=0.74, pad=0)

c=axlist[1].contourf(lon, lat, nosmokePM,60, transform=ccrs.PlateCarree())
axlist[1].set_title('nosmoke PM mort', fontsize=16)
cb1 = fig.colorbar(c, ax=axlist[1], orientation='horizontal', shrink=0.74, pad=0)

c=axlist[2].contourf(lon, lat, smokePMaf,60, transform=ccrs.PlateCarree())
axlist[2].set_title('smoke PM mort', fontsize=16)
cb = fig.colorbar(c, ax=axlist[2], orientation='horizontal', shrink=0.74, pad=0)

plt.savefig(out_fig_path + 'smoke_mort_cartopycheck_'+desc+'.png')
plt.show()
'''