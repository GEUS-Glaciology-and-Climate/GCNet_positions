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
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
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
df_all = pd.DataFrame()
i=0
# df_pos = df_pos.loc[df_pos.name_long=='JAR1',:]

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
    for d in pd.to_datetime(meta.loc[meta.ID==j, ['InstallationDate']].values[0]):
        s = tmp_interp.iloc[0:1,:]
        s.iloc[0,:2] = np.nan
        s.iloc[0,2]=[pd.to_datetime(str(d.year)+'-01-01')]
        tmp_interp = pd.concat((tmp_interp, s))
    for d in pd.to_datetime(meta.loc[meta.ID==j, ['LastValidDate']].values[0]):
        s = tmp_interp.iloc[0:1,:]
        s.iloc[0,:2] = np.nan
        s.iloc[0,2]=[pd.to_datetime(str(d.year)+'-12-31')]
        tmp_interp = pd.concat((tmp_interp, s))
    tmp_interp = tmp_interp.sort_values('date')
    tmp_interp = tmp_interp.loc[tmp_interp.date.notnull(), :].set_index('date')
    tmp_interp = tmp_interp.resample('H').mean()
    tmp_interp = tmp_interp.interpolate(method='spline',order=1, 
                                        limit_direction='both', fill_value="extrapolate")

    tmp_interp_y = tmp_interp.loc[[str(y)+'-06-01' for y in tmp_interp.index.year.unique()],:]

    # fig, ax = plt.subplots(2,1,sharex=True)
    # ax[0].plot(tmp_interp.index, tmp_interp.lon, marker='.', linestyle='None', label='hourly estimate')
    # ax[0].plot(tmp_interp_y.index, tmp_interp_y.lon, marker='d', linestyle='None', label='June 1st estimate')
    # ax[0].plot(tmp.index, tmp.lon, marker='o', linestyle='None', label='observed')
    # ax[0].set_ylabel('Longitude (deg W)')
    # ax[0].legend()
    # ax[0].grid()

    # ax[1].plot(tmp_interp.index, tmp_interp.lat, marker='.', linestyle='None')
    # ax[1].plot(tmp_interp_y.index, tmp_interp_y.lat, marker='d', linestyle='None')
    # ax[1].plot(tmp.index, tmp.lat, marker='o', linestyle='None')
    # ax[1].set_ylabel('Latitude (deg N)')
    # ax[1].grid()
    # plt.suptitle(station)
    
    tmp_interp.to_csv('output/'+station+'_position_interpolated.csv', float_format='%.5f')
                 
    tmp_interp_y['site'] = station
    if len(df_all)==0:
        df_all = tmp_interp_y[['site', 'lon','lat']].reset_index()
    else:
        df_all = pd.concat((df_all, tmp_interp_y[['site', 'lon','lat']].reset_index()))
    
    tmp_interp = tmp_interp.reset_index()
    tmp_interp_y = tmp_interp_y.reset_index()
    
    gdf = gpd.GeoDataFrame(tmp, geometry=gpd.points_from_xy(tmp.lon, tmp.lat))
    gdf = gdf.set_crs(4326)
    gdf = gdf.to_crs(3413)
    gdf['v'] = np.sqrt(gdf.geometry.x.diff()**2 + gdf.geometry.y.diff()**2) / (tmp.index.to_series().diff().dt.days/365)
    gdf['x'] = gdf.geometry.x/1000
    gdf['y'] = gdf.geometry.y/1000
    gdf2 = gpd.GeoDataFrame(tmp_interp_y, geometry=gpd.points_from_xy(tmp_interp_y.lon, tmp_interp_y.lat))
    gdf2 = gdf2.set_crs(4326)
    gdf2 = gdf2.to_crs(3413)
    gdf2['v'] = np.sqrt(gdf2.geometry.x.diff()**2 + gdf2.geometry.y.diff()**2) / (tmp_interp_y.date.diff().dt.days/365)
    gdf2['x'] = gdf2.geometry.x/1000
    gdf2['y'] = gdf2.geometry.y/1000
    print(j, station, '%0.0f'%gdf.loc[gdf.v.notnull(),'v'].median(), '%0.0f'%gdf2.loc[gdf2.v.notnull(),'v'].median())


    fig, ax=plt.subplots(1,1, figsize=(12,8))
    plt.subplots_adjust(left = 0.2, right=0.8)
    # ax=ax.flatten()
    ax = [ax]
    i = i+1
    ax[0].set_title(tmp.name_long.unique()[0],fontsize=14)
    plt.xlabel('Longitude ($^o$E)',fontsize=14)
    plt.ylabel('Latitude ($^o$N)',fontsize=14)
    w = tmp_interp_y.lon.max() - tmp.lon.min()
    xlim = [tmp_interp_y.lon.min()-w/6, tmp_interp_y.lon.max()+w/6]
    ax[0].set_xlim(xlim)
    h = tmp_interp_y.lat.max() - tmp_interp_y.lat.min()
    ylim = [tmp_interp_y.lat.min()-h/7, tmp_interp_y.lat.max()+h/8]
    ax[0].set_ylim(ylim)
    plt.plot(np.nan,np.nan, 'o',markerfacecolor='k', linestyle='None', label='GPS measurements')
    plt.plot(np.nan,np.nan, 'd',markerfacecolor='r',linestyle='None', label='inter- or extrapolated position on 1 June using \n spline fit on measured position')
    plt.ticklabel_format(axis='both',style='plain',useOffset=False)
    if station == 'Swiss Camp':
        loc = 'upper left'
    else:
        loc = 'best'
    plt.legend(loc=loc, fontsize=12)
    
    images = []
    for year in gdf.index.year.unique():
        tmp2 = tmp.loc[str(year),:].copy()
        ax[0].plot(tmp2.lon, tmp2.lat, 'k', markersize=10,
                   label='observed',
                   marker = 'o', linestyle='None')
        if tmp2.shape[0]>1:
            tmp2 = tmp2[['lat','lon']].resample('Y').mean()
        ax[0].annotate(str(tmp2.index.year.values[0]),
                       xy=(tmp2.lon, tmp2.lat),
                       xycoords='data',
                       xytext=(120, 0), 
                       textcoords='offset pixels',
                        fontsize=12,verticalalignment='center',
                        path_effects=[pe.withStroke(linewidth=4, foreground="white", alpha = 0.5)],
                        arrowprops=dict(arrowstyle="-", edgecolor='black'),
                        zorder=0)
        
        if make_gif == 1:
            filename='figs/'+df_pos.loc[df_pos.id==j,'name'].iloc[0]+'_'+str(year)+'.png'
            fig.savefig(filename, dpi=300)
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
                   xytext=(-250, 0), 
                   textcoords='offset pixels',
                   fontsize=11, color='gray',
                   arrowprops=dict(arrowstyle="-",edgecolor='lightgray'))
    scale = 500
    loc = 'lower right'
    if station in ['DYE-2', 'NASA-E', 'Tunu-N']:
        scale = 100
    if station in ['Humboldt']:
        scale = 100
    if station in ['NASA-SE']:
        scale = 10
        loc='lower left'
    label = str(scale)+' m'
    if station == 'Swiss Camp':
        scale = 1000
        label = '1 km'


    R = 6371e3 # metres
    phi = tmp_interp_y.lat.min() * np.pi/180 # φ, λ in radians
    d_lon = 1 * np.pi/180
    
    a = np.cos(phi) **2 * np.sin(d_lon/2) * np.sin(d_lon/2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    d = R * c # in km

    scalebar = AnchoredSizeBar(ax[0].transData,
                               scale/d, label, loc, 
                               pad=1,
                               color='k',
                               frameon=False,
                               size_vertical=(ylim[1]-ylim[0])/200,
                               fontproperties=fm.FontProperties(size=14),
                               )
    
    ax[0].add_artist(scalebar)
    filename='figs/'+df_pos.loc[df_pos.id==j,'name'].iloc[0]+'_final.png'
    fig.savefig(filename,dpi=300)
    if make_gif == 1:
        images.append(imageio.imread(filename))
        images.append(imageio.imread(filename))
        images.append(imageio.imread(filename))
        if station == 'Swiss Camp':
            imageio.mimsave('figs/gifs/'+station+'.gif', images,   duration=0.4)
        else:
            imageio.mimsave('figs/gifs/'+station+'.gif', images,   duration=0.6)

    df_all.to_csv('output/GC-Net_annual_summer_position_estimated.csv',index=None)


        