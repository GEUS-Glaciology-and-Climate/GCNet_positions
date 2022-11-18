# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
# import matplotlib.colormap as cm
import geopandas as gpd
from pandas.tseries.frequencies import to_offset
import matplotlib.patheffects as pe
try:
    url = "https://docs.google.com/spreadsheets/d/1R2SA7rqo9PHfAAGeSVgy7eWVHRugV8Z3nbWga5Xin1U/export?format=csv&gid=0"
    pd.read_csv(url).to_csv("data/GC-Net_yearly_positions.csv", index=None)
except:
    pass

df_pos = pd.read_csv("data/GC-Net_yearly_positions.csv")
plt.close('all')
# col = cm('Spectral',32)
abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'


i=0
for j in df_pos.id.unique():

    if df_pos.loc[df_pos.id==j,'name'].iloc[0] not in ['SWC', 'CP1', 'NAU', 
                                                       'GIT', 'TUN', 'DY2', 
                                                       'JR1', 'NAE','NSE', 'PET']:
        continue
    if j in [25, 26, 27, 28]:
        continue
    fig, ax=plt.subplots(1,1, figsize=(7,6))
    # ax=ax.flatten()
    ax = [ax]
    i = i+1
    tmp = df_pos.loc[df_pos.id==j,:].reset_index(drop=True)

    ax[0].set_title(abc[i-1]+'. '+tmp.name_long.unique()[0])

    tmp_interp= pd.DataFrame()
    tmp_interp['date'] = pd.to_datetime(np.append( tmp.date.values, ['1995-05-01', '2022-05-01']), errors='coerce')
    tmp_interp['lon'] = tmp.lon
    tmp_interp['lat'] = tmp.lat
    tmp_interp = tmp_interp.sort_values('date')
    tmp_interp = tmp_interp.loc[tmp_interp.date.notnull(), :].set_index('date')
    tmp_interp = tmp_interp.resample('3M').mean()
    tmp_interp=tmp_interp.interpolate(method='spline', order=1, limit_direction='both', fill_value="extrapolate")
    tmp_interp_y = tmp_interp.resample('Y').first()
    tmp_interp_y.index = tmp_interp_y.index - to_offset("11M")
    tmp_interp = tmp_interp.reset_index()
    tmp_interp_y = tmp_interp_y.reset_index()

    # ax[0].plot(tmp_interp.lon, tmp_interp.lat, marker = '.', color='tab:red', linestyle='None')
    ax[0].plot(tmp_interp_y.lon, tmp_interp_y.lat, 
               label='annual inter- or extrapolated',
               marker = 'd',color='tab:red',  linestyle='None')
    
    for k in tmp_interp_y.index:
        ax[0].text(tmp_interp_y.lon[k]+tmp_interp_y.lon.diff().abs().mean()*2,
                tmp_interp_y.lat[k], 
                tmp_interp_y.date[k].year, fontsize=11, color='lightgray')
    ax[0].plot(tmp.lon, tmp.lat, 
               label='observed',
               marker = 'o', linestyle='None')

    for k in tmp.index:
        ax[0].text(tmp.lon[k]+tmp_interp_y.lon.diff().abs().mean()*2,
                tmp.lat[k], 
                tmp.date[k], fontsize=16,
                path_effects=[pe.withStroke(linewidth=4, foreground="white")])
    xlim = ax[0].get_xlim()
    ax[0].set_xlim(xlim[0], xlim[1] + 0.2*(xlim[1] - xlim[0]))
    plt.legend()
    plt.xlabel('Longitude ($^o$E)')
    plt.ylabel('Latitude ($^o$E)')
    fig.savefig('figs/'+df_pos.loc[df_pos.id==j,'name'].iloc[0]+'.png')
    print('![](https://github.com/GEUS-Glaciology-and-Climate/GCNet_positions/blob/main/figs/'+df_pos.loc[df_pos.id==j,'name'].iloc[0]+'.png)')
    

        