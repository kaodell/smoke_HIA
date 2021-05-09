#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pool_asth_OR.py
    python script to plot asthma odds ratios, pool ORs across studies,
    and plot with those from Borchers-Arrigada meta-analysis
Created on Thu Sep 24 16:02:11 2020
@author: kodell
"""
#%% user inputs
# file locations
epi_rr_2pool_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/acute_analysis_info/WFS_Asth_RRs_2pool_final.csv'
# out file and figure path and descrption
out_fp = '/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA data/outputs/acute_outcomes/'
out_fig_path ='/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/acute_outcomes/'
out_desc = '_wMCweight_metamethod'

#%% load modules
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import pandas as pd
pio.renderers.default = "chrome"
pio.templates.default = "seaborn"

#%% load data
# asthma odds ratios to include in the pooled estimate
asthRR_df_2pool = pd.read_csv(epi_rr_2pool_fn)


#%% calculate standard error (SE) for each study
# have to convert to betas first
asthRR_df_2pool['beta'] = np.log(asthRR_df_2pool['OR 10ug/m^3'])/10.0
asthRR_df_2pool['beta_se'] = np.mean([(np.log(asthRR_df_2pool['OR 95% UCI 10ug/m^3'])/10.0)-asthRR_df_2pool['beta'],
                                        asthRR_df_2pool['beta']-(np.log(asthRR_df_2pool['OR 95% LCI 10ug/m^3'])/10.0)],axis=0)/1.96

#%% monte carlo pooling
# create array to save pooled ORs
fill_array = [-999,-999]
pool_OR = pd.DataFrame(data = {'Study':['pooled estimate','pooled estimate'],
                               'Region':['US','US'],
                               'Outcome(s)':['asth ED', 'asth Hosp'],
                               'OR 10ug/m^3':fill_array,
                               'OR 95% LCI 10ug/m^3':fill_array,
                               'OR 95% UCI 10ug/m^3':fill_array,
                               })

# pool using monte carlo
n_samples = 1000 # actual n samples is n_samples*1/se
# group by outcome
asthRR_df_2pool_outcome_groups = asthRR_df_2pool.groupby(by='Outcome(s)')
j=0
for outcome in ['asth ED', 'asth Hosp']:
    outcome_group = asthRR_df_2pool_outcome_groups.get_group(outcome)
    all_outcome_samples = np.array([])

    for i in range(outcome_group.shape[0]):
        mean = outcome_group['beta'].iloc[i]
        sd = outcome_group['beta_se'].iloc[i] # have to use standard error, don't have standard deviation or n
        w_samples = 1.0/sd
        samples = np.random.normal(loc=mean,scale=sd,size=int(n_samples*w_samples))
        all_outcome_samples = np.append(all_outcome_samples, samples)
        print(sd,w_samples)
    # calculate mean, 2.5 and 97.5 percentile for outcome samples
    outcome_ORpmc = np.exp(10.0*np.mean(all_outcome_samples))
    outcome_OR_pmc_lci = np.exp(10.0*np.percentile(all_outcome_samples,2.5))
    outcome_OR_pmc_uci = np.exp(10.0*np.percentile(all_outcome_samples,97.5))

    # add to array
    pool_OR['OR 10ug/m^3'].iloc[j] = outcome_ORpmc
    pool_OR['OR 95% LCI 10ug/m^3'].iloc[j] = outcome_OR_pmc_lci
    pool_OR['OR 95% UCI 10ug/m^3'].iloc[j] = outcome_OR_pmc_uci
    j+=1


# save the pooled values to a new spreadsheet
pool_OR.to_csv(out_fp + 'asth_COPD_ORs_wpool_'+out_desc+'.csv')

#%% make forrest plot with this data and that from Borchers Arriagada
colors_all = ['cornflowerblue','indianred']
colors_pool = ['darkblue','darkred']

# add columns of the error for plotting
asthRR_df_2pool['error_y'] = asthRR_df_2pool['OR 95% UCI 10ug/m^3'].values - asthRR_df_2pool['OR 10ug/m^3'].values
asthRR_df_2pool['error_y_minus'] = asthRR_df_2pool['OR 10ug/m^3'].values - asthRR_df_2pool['OR 95% LCI 10ug/m^3'].values

pool_OR['error_y'] = pool_OR['OR 95% UCI 10ug/m^3'].values - pool_OR['OR 10ug/m^3'].values
pool_OR['error_y_minus'] = pool_OR['OR 10ug/m^3'].values - pool_OR['OR 95% LCI 10ug/m^3'].values


fig = make_subplots(rows=1,cols=2,subplot_titles= ['Asthma Hospitalizaiton ORs',
                                                  'Asthma ED Visits ORs'],
                    x_title='Odds Ratio per 10 	&#956;g m<sup>-3</sup> PM<sub>2.5<sub>')
i = 0
col=1
for outcome in ['asth Hosp','asth ED']:
    # add meta-analysis from Borchers Arriagada et al. 2019
    if outcome == 'asth Hosp':
            fig.add_trace(go.Scatter(x = [1.08], y = ['Borchers Arriagada et al. (2019)'],
                             error_x=dict(symmetric=False,type='data',
                                          array = [1.14-1.08],
                                          arrayminus = [1.08-1.03]),
                      mode='markers',marker_color=colors_pool[col-1],marker_symbol='square',
                      marker_size=12),row=1,col=col)
    if outcome == 'asth ED':
            fig.add_trace(go.Scatter(x = [1.07], y = ['Borchers Arriagada et al. (2019)'],
                             error_x=dict(symmetric=False,type='data',
                                          array = [1.11-1.07],
                                          arrayminus = [1.07-1.03]),
                      mode='markers',marker_color=colors_pool[col-1],marker_symbol='square',
                      marker_size=12),row=1,col=col)

    # add pooled RR
    pool_ind = np.where(pool_OR['Outcome(s)']==outcome)[0]
    fig.add_trace(go.Scatter(x = pool_OR['OR 10ug/m^3'].iloc[pool_ind],
                             y = pool_OR['Study'].iloc[pool_ind],
                             error_x=dict(symmetric=False,type='data',
                                          array = pool_OR['error_y'].iloc[pool_ind],
                                          arrayminus = pool_OR['error_y_minus'].iloc[pool_ind]),
                      mode='markers',marker_color=colors_pool[col-1],
                      marker_size=12),row=1,col=col)

    # add odds ratios going into pooled OR
    ors2plot = asthRR_df_2pool_outcome_groups.get_group(outcome)
    fig.add_trace(go.Scatter(x = ors2plot['OR 10ug/m^3'], y = ors2plot['Study'],
                             error_x=dict(symmetric=False,type='data',
                                          array = ors2plot['error_y'],
                                          arrayminus = ors2plot['error_y_minus']),
                      mode='markers',marker_color=colors_all[col-1],
                      marker_size=12),row=1,col=col)
                             
    # add one line
    fig.add_trace(go.Scatter(x = [1]*(1+len(ors2plot['Study'])), 
                             y = np.append(ors2plot['Study'].values,
                                           ['Borchers Arriagada et al. (2019)']),
                      mode='lines',marker_color='grey'),row=1,col=col)
    col += 1
    i += 1
fig.update_traces(marker=dict(line=dict(width=2,)),selector=dict(mode='markers'))
fig.update_layout(plot_bgcolor = 'white',showlegend=False,font_family='Arial',font_size=16)
fig.update_xaxes(ticks="outside",tickvals= [0.95,1,1.05,1.10,1.15,1.20,1.25,1.30,1.35],
                 tickwidth=1, tickcolor='grey', ticklen=10,
                 showline=True, linewidth=1, linecolor='black')

fig.update_yaxes(title_text='Study',
                 ticks="outside", tickwidth=1, tickcolor='grey', ticklen=10,
                 showline=False, linewidth=1, linecolor='black',row=1,col=1)

fig.show()
fig.write_image('/Users/kodell/Local Google Drive /CSU/Research/NSF/smoke-specific HIA/smoke-specific HIA figures/acute_outcomes/forestplot_wborchers'+out_desc+'.png',
                width=1600,height=600)
