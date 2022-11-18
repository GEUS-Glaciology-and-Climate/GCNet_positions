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
for j in [1,2,3,4,5,6,7,8,9,10, 11,12,15,22,23,24]:

    if df_pos.loc[df_pos.id==j,'name'].iloc[0] in ['EGP', 'NEM', 'SDO', 'SDL', 'SUM', 'GITS']:
       continue   
    fig, ax=plt.subplots(1,3, figsize=(15,9))
    plt.subplots_adjust(left=0.08, right=0.99,
                        bottom=0.08,top=0.95,
                        hspace = 0.35, wspace=0.25)
    # ax=ax.flatten()
    # ax = [ax]
    i = i+1
    tmp = df_pos.loc[df_pos.id==j,:].reset_index(drop=True)

    ax[0].plot(tmp.lon, tmp.lat, marker = 'o', linestyle='None')
    
    texts = []
    for k in tmp.index:
        texts.append(ax[0].text(tmp.lon[k],
                                tmp.lat[k], 
                                tmp.date[k], fontsize=16,
                                path_effects=[pe.withStroke(linewidth=4, foreground="white")]))
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

    ax[0].plot(tmp_interp.lon, tmp_interp.lat, marker = '.', color='tab:red', linestyle='None')
    ax[0].plot(tmp_interp_y.lon, tmp_interp_y.lat, marker = 'd',color='tab:red',  linestyle='None')
    for k in tmp_interp_y.index:
        texts.append(ax[0].text(tmp_interp_y.lon[k],
                                tmp_interp_y.lat[k], 
                                tmp_interp_y.date[k].year, fontsize=12))
    # interpolating the x
    ax[1].set_title('latitude vs time')
    ax[1].plot(pd.to_datetime(tmp.date, errors='coerce'), tmp.lat, marker = 'o', linestyle='None')
    ax[1].plot(tmp_interp.date, tmp_interp.lat, marker = '.', color='tab:red', linestyle='None')
    ax[1].plot(tmp_interp_y.date, tmp_interp_y.lat, marker = 'd',color='tab:red',  linestyle='None')
    
    # interpolating the y
    ax[2].set_title('longitude vs time')
    ax[2].plot(pd.to_datetime(tmp.date, errors='coerce'), tmp.lon, marker = 'o', linestyle='None')
    ax[2].plot(tmp_interp.date, tmp_interp.lon, marker = '.', color='tab:red', linestyle='None')
    ax[2].plot(tmp_interp_y.date, tmp_interp_y.lon, marker = 'd',color='tab:red',  linestyle='None')


        