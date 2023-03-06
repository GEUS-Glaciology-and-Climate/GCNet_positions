# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from adjustText import adjust_text
# import matplotlib.colormap as cm
import geopandas as gpd
from pandas.tseries.frequencies import to_offset
import matplotlib.patheffects as pe
from datetime import datetime as dt
import time
import imageio.v2 as imageio

# (down)loading coordinates spreadsheet 
try:
    url = "https://docs.google.com/spreadsheets/d/1R2SA7rqo9PHfAAGeSVgy7eWVHRugV8Z3nbWga5Xin1U/export?format=csv&gid=0"
    pd.read_csv(url).to_csv("data/GC-Net_yearly_positions.csv", index=None)
except:
    print('Cannot access online file, using local file')
    pass
df_pos = pd.read_csv("data/GC-Net_yearly_positions.csv")
meta = pd.read_csv('data/GC-Net_location.csv', skipinitialspace=True)

make_gif = 0

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

plt.close('all')
# col = cm('Spectral',32)
abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

i=0
for j in  df_pos.id.unique():
    station=df_pos.loc[df_pos.id==j,'name_long'].iloc[0]
    if station in ['GITS', 'South Dome', 'Saddle', 'Summit','NEEM']:
        continue
    tmp = df_pos.loc[df_pos.id==j,:].reset_index(drop=True).copy()
    tmp.date = pd.to_datetime(tmp.date, errors='coerce')
    tmp = tmp.set_index('date')
    tmp = tmp.loc[tmp.index.notnull(),: ]
    if len(tmp.index.year.unique())<3:
        continue
    
    tmp_interp= pd.DataFrame()
    tmp_interp['lon'] = tmp.lon.values
    tmp_interp['lat'] = tmp.lat.values
    tmp_interp['date'] = tmp.index.values
    for d in pd.to_datetime(meta.loc[meta.ID==j, ['InstallationDate','LastValidDate']].values[0]):
        s = tmp_interp.iloc[0:1,:]
        s.iloc[0,:2] = np.nan
        s.iloc[0,2]=[d]
        tmp_interp = pd.concat((tmp_interp, s))
    tmp_interp = tmp_interp.sort_values('date')
    tmp_interp = tmp_interp.loc[tmp_interp.date.notnull(), :].set_index('date')
    tmp_interp = tmp_interp.resample('3M').mean()
    tmp_interp = tmp_interp.interpolate(method='spline', order=1, 
                                        limit_direction='both', fill_value="extrapolate")
    tmp_interp_y = tmp_interp.resample('Y').mean()
    tmp_interp_y.index = tmp_interp_y.index 
    tmp_interp = tmp_interp.reset_index()
    tmp_interp_y = tmp_interp_y.reset_index()
    
    gdf = gpd.GeoDataFrame(tmp, geometry=gpd.points_from_xy(tmp.lon, tmp.lat))
    gdf = gdf.set_crs(4326)
    gdf = gdf.to_crs(3413)
    gdf['v'] = np.sqrt(gdf.geometry.x.diff()**2 + gdf.geometry.y.diff()**2) / (tmp.index.to_series().diff().dt.days/365)
    
    gdf2 = gpd.GeoDataFrame(tmp_interp_y, geometry=gpd.points_from_xy(tmp_interp_y.lon, tmp_interp_y.lat))
    gdf2 = gdf2.set_crs(4326)
    gdf2 = gdf2.to_crs(3413)
    gdf2['v'] = np.sqrt(gdf2.geometry.x.diff()**2 + gdf2.geometry.y.diff()**2) / (tmp_interp_y.date.diff().dt.days/365)
    print(j, station, '%0.0f'%gdf.loc[gdf.v.notnull(),'v'].median(), '%0.0f'%gdf2.loc[gdf2.v.notnull(),'v'].median())
    
    fig, ax=plt.subplots(1,1, figsize=(7,6))
    plt.subplots_adjust(right=0.8)
    # ax=ax.flatten()
    ax = [ax]
    i = i+1
    ax[0].set_title(tmp.name_long.unique()[0],fontsize=14)
    plt.xlabel('Longitude ($^o$E)',fontsize=14)
    plt.ylabel('Latitude ($^o$N)',fontsize=14)
    w = tmp_interp_y.lon.max() - tmp.lon.min()
    xlim = [tmp_interp_y.lon.min()-w/8, tmp_interp_y.lon.max()+w/8]
    ax[0].set_xlim(xlim)
    h = tmp_interp_y.lat.max() - tmp_interp_y.lat.min()
    ylim = [tmp_interp_y.lat.min()-h/8, tmp_interp_y.lat.max()+h/8]
    ax[0].set_ylim(ylim)
    plt.plot(np.nan,np.nan, 'o',markerfacecolor='k', linestyle='None', label='GPS measurements')
    plt.plot(np.nan,np.nan, 'd',markerfacecolor='r',linestyle='None', label='annual inter- or extrapolation')
    plt.legend(loc='upper left',fontsize=12)
    
    images = []
    for year in tmp.index.year.unique():
        tmp2 = tmp.loc[str(year),:].copy()
        ax[0].plot(tmp2.lon, tmp2.lat, 'k', markersize=10,
                   label='observed',
                   marker = 'o', linestyle='None')
        if tmp2.shape[0]>1:
            tmp2 = tmp2[['lat','lon']].resample('Y').mean()
        ax[0].annotate(str(tmp2.index.year.values[0]),
                       xy=(tmp2.lon, tmp2.lat),
                       xycoords='data',
                       xytext=(60, 0), 
                       textcoords='offset pixels',
                        fontsize=12,verticalalignment='center',
                        path_effects=[pe.withStroke(linewidth=4, foreground="white", alpha = 0.5)],
                        arrowprops=dict(arrowstyle="-", edgecolor='black'),
                        zorder=0)
        if make_gif == 1:
            filename='figs/'+df_pos.loc[df_pos.id==j,'name'].iloc[0]+'_'+str(year)+'.png'
            fig.savefig(filename)
            images.append(imageio.imread(filename))
            os.remove(filename)
    ax[0].plot(tmp_interp_y.lon, tmp_interp_y.lat, 
               label='annual inter- or extrapolated',
               marker = 'd',color='tab:red',  linestyle='None',
               zorder=0)
    
    for k in range(tmp_interp_y.shape[0]):
        ax[0].annotate(tmp_interp_y.date[k].year,
                   xy=(tmp_interp_y.lon[k], tmp_interp_y.lat[k]),
                   xycoords='data',
                   xytext=(-60, 0), 
                   textcoords='offset pixels',
                   fontsize=11, color='gray',
                   arrowprops=dict(arrowstyle="-",edgecolor='lightgray'))

    filename='figs/'+df_pos.loc[df_pos.id==j,'name'].iloc[0]+'_final.png'
    fig.savefig(filename)
    if make_gif == 1:
        images.append(imageio.imread(filename))
        images.append(imageio.imread(filename))
        images.append(imageio.imread(filename))
        if station == 'Swiss Camp':
            imageio.mimsave(station+'.gif', images,   duration=0.4)
        else:
            imageio.mimsave(station+'.gif', images,   duration=0.6)
    

        