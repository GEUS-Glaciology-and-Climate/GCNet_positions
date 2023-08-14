#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: jeb@geus.dk

obtain monthly coordinates from transmissions

input:
    aws IMEI numbered decoded SBD transmissions
    
output:
    monthly Google Earth placemarkers
    the attached table, equivalent with that for PROMICE ESSD positions table after https://github.com/GEUS-Glaciology-and-Climate/PROMICE_positions
    graphics like in the ESSD

"""
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib.pyplot import figure
import matplotlib.dates as mdates
from datetime import datetime
import sys
import simplekml
# import geopy.distance

# -------------------------------- set the working path automatically
if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/AWS/GCNet/GCNet_positions/'
os.chdir(base_path)
sys.path.append(base_path)

from datetime import datetime
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from datetime import date
import geopandas as gpd
from pyproj import Proj, transform
from PIL import Image
# from datetime import date
from datetime import timedelta
import ftplib
import calendar


# ----------------------------------------------------------
# ----------------------------------------------------------
# some main user-defined variables
plot_individual=1 ; site='JAR_O' #'SWC_O'# 'QAS_U' 
do_plot=1 # set to 1 if user wants plots
plt_map=0 # if set to 1 draws a map to the right of the graphic
ly='p' # either 'x' or 'p', 'x' is for display to local plots window, 'p' writes a .png graphic
do_NRT=1 # do near realtime? 1 f or yes
do_ftp=0 # set to 1 if push values
plot_stars_on_extreme_values=1 # like it says ;-)
open_fig_testing=0
n_std=1.96 # 1.96 sigma corresponds to 95%ile
min_years_of_record_for_robust_stats=3
plt_last_val_text=1

# get today's date used in naming graphics and help with NRT
today = date.today()
today = pd.to_datetime(date.today())#-timedelta(days=1) 
versionx= today.strftime('%Y-%m-%d')
current_month=today.strftime('%m')
current_year=today.strftime('%Y')

# graphics definitions
th=2 # line thickness
formatx='{x:,.3f}' ; fs=18
plt.rcParams["font.size"] = fs
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.grid'] = False
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.5
plt.rcParams['grid.color'] = "#C6C6C6"
plt.rcParams["legend.facecolor"] ='w'
plt.rcParams["mathtext.default"]='regular'
plt.rcParams['grid.linewidth'] = th/2
plt.rcParams['axes.linewidth'] = 1

# ## change to your system's login name to change dir for local work
if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/'
os.chdir(base_path)

meta = pd.read_csv('/Users/jason/Dropbox/AWS/GEUS_AWS_NRT_alerts/ancil/AWS_station_locations.csv')
meta = meta.rename({'stid': 'name'}, axis=1)
# meta.drop(meta[meta.name=='Roof_GEUS'].index, inplace=True) # drop original time column
# meta.drop(meta[meta.name=='Roof_PROMICE'].index, inplace=True) # drop original time column
# meta.drop(meta[meta.name=='UWN'].index, inplace=True) # drop original time column
# meta.drop(meta[meta.name=='NUK_U'].index, inplace=True) # drop original time column
# drop some sites from the list, sites not transmitting in the interest of time
# names=['HUM','DY2','CEN1','CEN2','CP1','NAE','NAU','NEM','NSE','SDM','SDL']
names=['HUM','JAR','NAU','SUM','CEN1','CEN2','QAS_Uv3','SWC','QAS_Uv3','NUK_U','UWN','Roof_PROMICE','Roof_GEUS','LYN_T','LYN_L','KPC_Lv3','KPC_Uv3','THU_L2','WEG_B','ZAK_Uv3','MIT'] #
for name in names:
    meta.drop(meta[meta.name==name].index, inplace=True) # drop original time column

meta=meta.sort_values(by='name', ascending=True)
# print(meta.columns)

names=meta.name

timeframe='day'
# timeframe='hour'
iyear=2021 ; fyear=2023
n_years=fyear-iyear+1

# years=np.arange(iyear,fyear+1).astype(int)

# today = date.today()
last_day= today.strftime('%d')

if do_NRT:
    target_year=2023
    t0=datetime(1900,1,1)
    t1=datetime(1900,int(current_month),int(last_day))
else:
    # target_year=2022
    # t0=datetime(1900,8,28)
    # t1=datetime(1900,10,2)
    target_year=2022
    t0=datetime(1900,8,15)
    t1=datetime(1900,3,2) #DY2 test
    t1=datetime(1900,2,27) #CP1 test
    t1=datetime(1900,10,15) 

n_days=365
if calendar.isleap(target_year):n_days=366
print(n_days)

anomaly_message=[] # empty list

# loop over all sites by uncommentingthe line "if i>=0:" and commenting out all others
# or choose a site for testing, for example "if name=='NUK_Uv3':"
n_sites=len(names)

if plot_individual:
    names=names[names==site]
    
for i,name in enumerate(names):
    if i>=0:

        print()
        print(i,n_sites-i,name)#[i],meta.name[i],names[i][0:5])
        
        site=name

        # fn=aws_data_path+site+'/'+site+'_day.csv'
        # df=pd.read_csv(fn)        
        # print(site)
        url = "https://thredds.geus.dk/thredds/fileServer/aws_l3_station_csv/level_3/{}/{}_{}.csv".format(site,site,timeframe)
        df = pd.read_csv(url)
        print(df)
        
        print(df.columns)

        # df.t_u
        # position
        n=20
        lat=np.nanmean(df.gps_lat)#[-n:])
        lon=-abs(np.nanmean(df.gps_lon))#[-n:])
        elev=np.nanmean(df.gps_alt)
        # print(lat,lon,elev)
        
        # df.index = pd.to_datetime(df.time)
        # fig, ax = plt.subplots(figsize=(10,10))
        # plt.plot(df.batt_v)
        # plt.ylabel('VDC')
        # plt.title(name+' 2023 from Thredds')
        # plt.setp(ax.xaxis.get_majorticklabels(), rotation=90,ha='center' )
#%%
        df["date"] = pd.to_datetime(df.time)

        df['year'] = pd.DatetimeIndex(df["date"]).year
        df['day'] = pd.DatetimeIndex(df["date"]).day
        df['hour'] = pd.DatetimeIndex(df["date"]).hour
        df['month'] = pd.DatetimeIndex(df["date"]).month
        df['doy'] = pd.DatetimeIndex(df["date"]).dayofyear
        df['doy_dec'] = df['doy']+(df['hour']-1)/23
        df.index = pd.to_datetime(df.time)
    
        # print(df.columns)
        # print(df)
        

    
        #----- compute time dependence of position
        
        i_year=2021
        f_year=2023
        n_years=f_year-i_year+1
        
        means=np.zeros(((4,n_years,12)))*np.nan
        counts=np.zeros(((3,n_years,12)))
        years=np.zeros(n_years*12)
        dec_year=np.zeros(n_years*12)
        months=np.zeros(n_years*12)
        lat_mean=np.zeros(n_years*12)
        lat_count=np.zeros(n_years*12)
        lon_mean=np.zeros(n_years*12)
        lon_count=np.zeros(n_years*12)
        elev_mean=np.zeros(n_years*12)
        elev_std=np.zeros(n_years*12)
        elev_count=np.zeros(n_years*12)
    
        
        cc=0
        for yy in range(n_years):
            print(yy+i_year)
            for mm in range(12):
                v=((df.year==yy+i_year)&(df.month==mm+1)&(np.isfinite(df['gps_lat'])))
                means[0,yy,mm]=np.nanmean(df['gps_lat'][v])
                counts[0,yy,mm]=np.sum(v)
                v=((df.year==yy+i_year)&(df.month==mm+1)&(np.isfinite(df['gps_lon'])))
                means[1,yy,mm]=-np.nanmean(df['gps_lon'][v])
                counts[1,yy,mm]=np.sum(v)
                v=((df.year==yy+i_year)&(df.month==mm+1)&(np.isfinite(df['gps_alt'])))
                means[2,yy,mm]=np.nanmean(df['gps_alt'][v])
                counts[2,yy,mm]=np.sum(v)
                means[3,yy,mm]=np.nanstd(df['gps_alt'][v])
                print(yy+i_year,mm,means[0,yy,mm],means[1,yy,mm],means[2,yy,mm])
                
    
                years[cc]=yy+i_year
                months[cc]=mm+1
                dec_year[cc]=yy+i_year+(mm/12)
                lat_count[cc]=counts[0,yy,mm]
                lon_count[cc]=counts[1,yy,mm]
                elev_count[cc]=counts[2,yy,mm]
                thresh=0.5
                n_max_poss_data=31.5
                if lat_count[cc]>n_max_poss_data*thresh:
                    lat_mean[cc]=means[0,yy,mm]
                    kml = simplekml.Kml(open=1)
                    namex=site+'_'+str(yy+i_year)+'_'+str(mm+1).zfill(2)
                    pnt = kml.newpoint(name=namex)
                    pnt.coords=[(-means[1,yy,mm],means[0,yy,mm])]
                    print(-means[1,yy,mm],means[0,yy,mm])
                    pnt.style.labelstyle.color = simplekml.Color.red  # Make the text red
                    pnt.style.labelstyle.scale = 1.5  # text scaling multiplier
                    pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
                    pnt.altitudemode = simplekml.AltitudeMode.relativetoground
                    
                    kml_ofile='./output/kml/monthly/'+namex+".kml"
                    kml.save(kml_ofile)
    
                if lon_count[cc]>n_max_poss_data*thresh:lon_mean[cc]=means[1,yy,mm]
                if elev_count[cc]>n_max_poss_data*thresh:elev_mean[cc]=means[2,yy,mm]
                if elev_count[cc]>n_max_poss_data*thresh:elev_std[cc]=means[3,yy,mm]
     
        
                cc+=1
        #%%
        lat_mean[lat_mean==0]=np.nan
        lon_mean[lon_mean==0]=np.nan
        
        
        # data=[[years],[means],[stds]]
        df1d = pd.DataFrame(columns = ['year', 'month','lat','lat count','lon','lon count','elev','elev std','elev count'])
        df1d.index.name = 'index'
        df1d["year"]=pd.Series(years)
        df1d['year'] = df1d['year'].apply(lambda x: '%.0f' % x)
        df1d["month"]=pd.Series(months)
        df1d['month'] = df1d['month'].apply(lambda x: '%.0f' % x)
        df1d["lat"]=pd.Series(lat_mean)
        df1d['lat'] = df1d['lat'].apply(lambda x: '%.5f' % x)
        df1d["lat count"]=pd.Series(lat_count)
        df1d["lon"]=pd.Series(-lon_mean)
        df1d['lon'] = df1d['lon'].apply(lambda x: '%.6f' % x)
        df1d["lon count"]=pd.Series(lon_count)
        df1d["elev"]=pd.Series(elev_mean)
        df1d['elev'] = df1d['elev'].apply(lambda x: '%.2f' % x)
        df1d["elev std"]=pd.Series(elev_std)
        df1d['elev std'] = df1d['elev std'].apply(lambda x: '%.2f' % x)
        df1d["elev count"]=pd.Series(elev_count)
        ofile='./output/'+site+'_positions_monthly'
        df1d.to_csv(ofile+'.csv',index=None)
        print(df1d)
    # df1d.to_excel(ofile+'.xlsx')
    
    # v=df.Latitude[np.isfinite(df.Latitude)]
    # lat=np.nanmean(df.gps_lat[:-24])
    # lon=-np.nanmean(df.gps_lon[:-24])

    # elevx=np.nanmean(df.gps_alt[:-24])
    
    # print(lat,lon,elev)

        

#%%

# # graphics layout option
# ly='p' # p for .png, x for console only
# wo=1 # write out stats
# write_fig=1 # write out figures

# sites=['SWC','SDM','NSE','CP1','NAU','CEN','NEM','JAR']
# # sites=['SWC']

# n_sites=len(sites)

# stats=np.zeros((7,n_sites))
# date0=['']*n_sites
# date1=['']*n_sites
# years=np.zeros(n_sites)

# for site_index,site in enumerate(sites):
#     if site_index!=10:
# #     # if site=='SWC':
#     # if site=='CP1':
#     # if site=='SDM':

#         df=pd.read_csv('./output/'+site+'_positions_monthly.csv')
#         for col in df.columns:
#             df[col] = pd.to_numeric(df[col])
#         df['day']=15
#         df['date']=pd.to_datetime(df[["year", "month","day"]])
        
#         df['elev'][df['elev']==0]=np.nan

#         v=np.where(~np.isnan(df['lat']))
#         v=v[0]
#         n=len(v)
#         print(v[0],df['date'][v[0]])
#         x=[df['date'][v[0]],df['date'][v[n-1]]]

#         lat=df.lat
#         lon=df.lon
#         elev=df.elev
#         st=site_index

#         plt.close()
#         fig, ax = plt.subplots(3,1,figsize=(10,15))


#         cc=0
#         stnam=site
#         if ly!='n':
#             ax[cc].plot(df['date'],df['lat'],'.',color='k')
#             ax[0].set_title(stnam+' latitude')
#             ax[cc].set_xticklabels([])

#             cc+=1
#             ax[cc].plot(df['date'],df['lon'],'.',color='k')
#             ax[cc].set_title(stnam+' longitude')
#             ax[cc].set_xticklabels([])

#             cc+=1
#             ax[cc].plot(df['date'],df['elev'],'.',color='k')
#             ax[cc].set_ylabel('meters a.s.l')
        

#         cc=0        
#         stats[0,st]=lat[v[0]] 
#         ax[cc].plot(df['date'][v[0]],stats[0,st],'s',color='g',alpha=0.5,label='first valid datum')        
#         stats[1,st]=lat[v[n-1]]
#         ax[cc].plot(df['date'][v[n-1]],stats[1,st],'s',color='r',alpha=0.5,label='last valid datum')
#         y=[stats[0,st],stats[1,st]]
#         ax[cc].plot(x, y, '--',c='b',label='approximation')
#         ax[cc].legend()

#         cc+=1

#         stats[2,st]=lon[v[0]] 
#         ax[cc].plot(df['date'][v[0]],stats[2,st],'s',color='g',alpha=0.5)        
#         stats[3,st]=lon[v[n-1]]
#         ax[cc].plot(df['date'][v[n-1]],stats[3,st],'s',color='r',alpha=0.5)

#         y=[stats[2,st],stats[3,st]]
#         ax[cc].plot(x, y, '--',c='b')

#         cc+=1

#         stats[4,st]=elev[v[0]] 
#         ax[cc].plot(df['date'][v[0]],stats[4,st],'s',color='g',alpha=0.5)        
#         stats[5,st]=elev[v[n-1]]
#         ax[cc].plot(df['date'][v[n-1]],stats[5,st],'s',color='r',alpha=0.5)
        
#         plt.setp(ax[cc].xaxis.get_majorticklabels(), rotation=90,ha='center' )
#         # ax[cc].xaxis.set_major_formatter(mdates.DateFormatter('%Y %b %d'))

#         de=stats[5,st]-stats[4,st]

#         ax[2].set_title(stnam+' elevation, change = '+str('%.0f'%de)+' m')

#         y=[stats[4,st],stats[5,st]]
#         ax[cc].plot(x, y, '--',c='b')
        
#         coords_1 = (stats[0,st],stats[2,st])
#         coords_2 = (stats[1,st],stats[3,st])
#         stats[6,st]=geopy.distance.distance(coords_1, coords_2).m
#         # print(df.name[k],df.name[k],str("%.0f"%(dist)))
        
#         date0[st]=df['date'][v[0]]
#         date1[st]=df['date'][v[n-1]]
        
#         duration = df['date'][v[n-1]] - df['date'][v[0]]                         # For build-in functions
#         duration_in_s = duration.total_seconds()
#         years[st] = divmod(duration_in_s, 86400)[0]/365.25
#         print(years[st])


#         if ((ly=='p')&(write_fig)): plt.savefig('./figs/'+stnam+'.png', dpi=72,bbox_inches = 'tight')
#         if ly=='x':plt.show()

# #%% output multi-site table

# df2 = pd.DataFrame(columns=['site name','first valid date','latest valid date','delta time','first valid latitude, °N','latest valid latitude, °N','first valid longitude, °W','latest valid longitude, °W','displacement, m','displacement rate, m/y','first valid elevation, m','latest valid elevation, m','elevation change, m'])
# df2['site name']=pd.Series(sites)
# df2['first valid date']=pd.Series(date0[:])
# df2['latest valid date']=pd.Series(date1[:])
# df2['first valid latitude, °N']=pd.Series(stats[0,:])
# df2['latest valid latitude, °N']=pd.Series(stats[1,:])
# df2['first valid longitude, °W']=pd.Series(stats[2,:])
# df2['latest valid longitude, °W']=pd.Series(stats[3,:])
# df2['first valid elevation, m']=pd.Series(stats[4,:])
# df2['latest valid elevation, m']=pd.Series(stats[5,:])
# df2["displacement, m"]=pd.Series(stats[6,:])
# df2['elevation change, m']=df2['latest valid elevation, m']-df2['first valid elevation, m']
# df2["delta time"]=pd.Series(years[:])

# df2["displacement rate, m/y"]=df2["displacement, m"]/df2["delta time"]
# df2["displacement rate, m/y"][df2['site name']=='KAN_B']=np.nan
# df2['displacement, m'] = df2['displacement, m'].map(lambda x: '%.0f' % x)
# df2['elevation change, m'] = df2['elevation change, m'].map(lambda x: '%.1f' % x)
# df2['delta time'] = df2['delta time'].map(lambda x: '%.1f' % x)
# df2['first valid elevation, m'] = df2['first valid elevation, m'].map(lambda x: '%.0f' % x)
# df2['latest valid elevation, m'] = df2['latest valid elevation, m'].map(lambda x: '%.0f' % x)
# df2['first valid latitude, °N'] = df2['first valid latitude, °N'].map(lambda x: '%.6f' % x)
# df2['latest valid latitude, °N'] = df2['latest valid latitude, °N'].map(lambda x: '%.6f' % x)
# df2['first valid longitude, °W'] = df2['first valid longitude, °W'].map(lambda x: '%.6f' % x)
# df2['latest valid longitude, °W'] = df2['latest valid longitude, °W'].map(lambda x: '%.6f' % x)
# df2['displacement rate, m/y'] = df2['displacement rate, m/y'].map(lambda x: '%.1f' % x)

# if wo:df2.to_csv('./output/GCN_positions_distance_stats.csv',sep=';')
# if wo:df2.to_excel('./output/GCN_positions_distance_stats.xlsx')

# print("average displacement rate",np.nanmean(df2["displacement rate, m/y"].astype(float)))
# print("average displacement rate",np.nanstd(df2["displacement rate, m/y"].astype(float)))
# print("elevation change, m",np.nanmean(df2["elevation change, m"].astype(float)))
# print("elevation change, m",np.nanstd(df2["elevation change, m"].astype(float)))

