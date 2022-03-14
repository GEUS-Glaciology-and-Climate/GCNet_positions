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
import geopy.distance

# -------------------------------- set the working path automatically
if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/AWS/GCNet_positions/'
os.chdir(base_path)
sys.path.append(base_path)

figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
th=1 
font_size=18
# plt.rcParams['font.sans-serif'] = ['Georgia']
plt.rcParams["font.size"] = font_size
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 1
plt.rcParams['grid.color'] = "#cccccc"
plt.rcParams["legend.facecolor"] ='w'
plt.rcParams["mathtext.default"]='regular'
plt.rcParams['grid.linewidth'] = th
plt.rcParams['axes.linewidth'] = th #set the value globally
plt.rcParams['figure.figsize'] = 17, 10
plt.rcParams["legend.framealpha"] = 0.8
plt.rcParams['figure.figsize'] = 5, 4

font_size=20

sites=['SWC','SDM','NSE','CP1','NAU','CEN','NEM','JAR']
# sites=['SWC']
# sites=['JAR']

n_sites=len(sites)

for site_index,site in enumerate(sites):
    if site_index!=10:
#     # if site=='SWC':
    # if site=='CP1':
    # if site=='SDM':
    # if site=='CEN':

        if site=='SWC':
            site2='Swiss Camp'
            IMEI='300534060524310'# ; lat=69+33.21479/60 ;lon=-(49+22.297377/60); elev=1148
            swd_coef=12.00 ; swu_coef=18.00
        
        if site=='NUK_U':
            site2='NUK_U'
            xtit='UTC time, year 2021'
            swd_coef=1 ; swu_coef=1
            
        if site=='SDM':
            IMEI='300534060529220' ; site2='South Dome' #; lat=63.14889	;lon=-44.81667; elev=2895
            swd_coef=14.75 ; swu_coef=14.44
        
        if site=='NSE':
            IMEI='300534062416350' ; site2='NASA-SE' #; lat=66.2866629	;lon=-42.2961828; elev=2386
            swd_coef=12.66 ; swu_coef=13.9
            swd_coef=12.66 ; swu_coef=13.90
        
        if site=='CP1':
            IMEI='300534062418810' ; site2='Crawford Pt. (CP1)'#; lat=69.8819; lon=-46.97358; elev=1958
            swd_coef=15.72 ; swu_coef=12.18
                
        if site=='CEN':
            IMEI='300534062419950'; site2='Camp Century' #; lat=69.8819; lon=-46.97358; elev=1958
        
        if site=='NAU':
            IMEI='300534062412950'; site2='NASA-U' #; lat=69.8819; lon=-46.97358; elev=1958

        if site=='NEM':
            IMEI='300534062417910' ; site2='NEEM' #; lat=77.5022		;lon=-50.8744; elev=2454

        if site=='JAR':
            IMEI='300534060520350' ;site2='JAR'
        
        # if site=='QAS_U':
        #     site2='QAS_U'
        #     xtit='UTC time, year 2021'
        print(site)
        # -------------- copy of raw file names to an easier to understand format
        if ((site!='NUK_U')&(site!='QAS_U')):
            fn='/Users/jason/Dropbox/AWS/aws_data/AWS_'+IMEI+'.txt' 
        else:
            fn='/Users/jason/Dropbox/AWS/PROMICE/PROMICE_v04/'+site+'_hour_v04.txt'
        # os.system('/bin/cp '+fn+' ./AWS_info/AWS_data/'+site+'.txt')
        
        # fn='./AWS_info/AWS_data/'+site+'.txt'
        # cols=np.arange(0,49).astype(str)
        
        cols=['time','counter','Pressure_L','Pressure_U','Asp_temp_L','Asp_temp_U','Humidity_L','Humidity_U','WindSpeed_L','WindDirection_L','WindSpeed_U','WindDirection_U','SW Downward','SW Upward','LW Downward','LW Upward','TemperatureRadSensor','SR_L','SR_U','T_firn_1','T_firn_2','T_firn_3','T_firn_4','T_firn_5','T_firn_6','T_firn_7','T_firn_8','T_firn_9','T_firn_10','T_firn_11','Roll','Pitch','Heading','Rain_amount_L','Rain_amount_U','counterx','Latitude','Longitude','Altitude','ss','Giodal','GeoUnit','Battery','NumberSatellites','HDOP','FanCurrent_L','FanCurrent_U','Quality','LoggerTemp']
        varnamx=['time','counter','Pressure_L','Pressure_U','air temperature','air temperature','Humidity_L','Humidity_U','WindSpeed_L','WindDirection_L','WindSpeed_U','WindDirection_U','SW Downward','SW Upward','LW Downward\nSky Tempearure_effective','LW Upward','TemperatureRadSensor','SR_L','SR_U','T_firn_1','T_firn_2','T_firn_3','T_firn_4','T_firn_5','T_firn_6','T_firn_7','T_firn_8','T_firn_9','T_firn_10','T_firn_11','Roll','Pitch','Heading','rainfall','rainfall','counterx','Latitude','Longitude','Altitude','ss','Giodal','GeoUnit','Battery','NumberSatellites','HDOP','FanCurrent_L','FanCurrent_U','Quality','LoggerTemp']
        unitsx=['time','counter','Pressure_L','Pressure_U','deg. C','deg. C','Humidity_L','Humidity_U','WindSpeed_L','WindDirection_L','WindSpeed_U','WindDirection_U','SW Downward','SW Upward','LW Downward','LW Upward','TemperatureRadSensor','SR_L','SR_U','T_firn_1','T_firn_2','T_firn_3','T_firn_4','T_firn_5','T_firn_6','T_firn_7','T_firn_8','T_firn_9','T_firn_10','T_firn_11','Roll','Pitch','Heading','mm','mm','counterx','Latitude','Longitude','Altitude','ss','Giodal','GeoUnit','Battery','NumberSatellites','HDOP','FanCurrent_L','FanCurrent_U','Quality','LoggerTemp']
        
        
        cols[0]='time'
        skip=1
        if site=='NSE':
            skip=3
            df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
        if site=='NEM':
            df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
        if site=='SDM':
            df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
        if site=='NAU':
            skip=4
            df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
        if site=='CP1':
            df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
        if site=='CEN':
            skip=4
            df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
        if site=='JAR':
            df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)

            # cols=['time','counter','Pressure_L','Pressure_U','Asp_temp_L','Asp_temp_U','Humidity_L','Humidity_U','WindSpeed_L','WindDirection_L','WindSpeed_U','WindDirection_U','SW Downward','SW Upward','LW Downward','LW Upward','TemperatureRadSensor','SR_L','SR_U','T_firn_1','T_firn_2','T_firn_3','T_firn_4','T_firn_5','T_firn_6','T_firn_7','T_firn_8','T_firn_9','T_firn_10','T_firn_11','Roll','Pitch','Heading','Rain_amount_L','Rain_amount_U','counterx','Latitude','Longitude','Altitude','ss','Giodal','GeoUnit','Battery','NumberSatellites','HDOP','FanCurrent_L','FanCurrent_U','Quality','LoggerTemp']
        #     df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
            
        if ((site=='NUK_U')or(site=='QAS_U')):
            # df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
            df=pd.read_csv(fn, delim_whitespace=True)
            df.columns
            df[df == -999.0] = np.nan
            # df[df=="NAN"]=np.nan
            # for col in df.columns:
            #     df[col] = pd.to_numeric(df[col])
        
            # df.rename({'AirPressure(hPa)': 'Pressure_L'}, axis=1, inplace=True)
            # df.rename({'RelativeHumidity(%)': 'Humidity_L'}, axis=1, inplace=True)
        
            # df.rename({'HeightStakes(m)': 'SR_L'}, axis=1, inplace=True)
            # df.rename({'HeightSensorBoom(m)': 'SR_U'}, axis=1, inplace=True)
            # df.rename({'AirTemperature(C)': 'Asp_temp_L'}, axis=1, inplace=True)
            # df.rename({'Ac_Rain(mm)': 'Rain_amount_L'}, axis=1, inplace=True)
        
            # df.rename({'LongwaveRadiationDown(W/m2)': 'LWD'}, axis=1, inplace=True)
            # df.rename({'LongwaveRadiationUp(W/m2)': 'LWU'}, axis=1, inplace=True)
            # df.rename({'ShortwaveRadiationDown_Cor(W/m2)': 'SW Downward'}, axis=1, inplace=True)
            # df.rename({'ShortwaveRadiationUp_Cor(W/m2)': 'SW Upward'}, axis=1, inplace=True)
            # df.rename({'ElevationGPS(m)': 'Altitude'}, axis=1, inplace=True)
            # df.rename({'LongitudeGPS(degW)': 'Longitude'}, axis=1, inplace=True)
            # df.rename({'LatitudeGPS(degN)': 'Latitude'}, axis=1, inplace=True)
            # df.rename({'WindSpeed(m/s)': 'WindSpeed_L'}, axis=1, inplace=True)
        
            df.rename({'Year': 'year'}, axis=1, inplace=True)
            df.rename({'MonthOfYear': 'month'}, axis=1, inplace=True)
            df.rename({'DayOfMonth': 'day'}, axis=1, inplace=True)
            df.rename({'HourOfDay(UTC)': 'hour'}, axis=1, inplace=True)
        
            df["time"]=pd.to_datetime(df[['year', 'month', 'day', 'hour']])
            df.index = pd.to_datetime(df.time)
        
            # print(df.columns)
            # print(df)
            
        if site=='SWC'or site=='JAR':
            skip=4
            cols=['time','seconds_since_1990','Pressure_L','Asp_temp_L','Humidity_L','WindSpeed_L','winddirection_s_l',
                  'SW Upward','SW Downward',
                 'LW Downward','LW Upward','TemperatureRadSensor','SR_L','SR_U','solar?',
                 'thermistorstring_1','thermistorstring_2','thermistorstring_3','thermistorstring_4','thermistorstring_5','thermistorstring_6','thermistorstring_7','thermistorstring_8',
                 'roll','pitch','heading','Rain_amount_L','gpstime','Latitude','Longitude','Altitude','giodal','geounit?',
                 'battvolt','?1','asp_temp_u','humidity_u','##','##2','##3']
            varnamx=['time','seconds_since_1990','pressure_l','Asp_temp_L','humidity_l','WindSpeed_L','winddirection_s_l','swupper','swlower',
                 'lwupper','lwlower','temperatureradsensor','sr_l','sr_u','solar?',
                 'thermistorstring_1','thermistorstring_2','thermistorstring_3','thermistorstring_4','thermistorstring_5','thermistorstring_6','thermistorstring_7','thermistorstring_8',
                 'roll','pitch','heading','Rain_amount_L','gpstime','Latitude','Longitude','Altitude','giodal','geounit?',
                 'battvolt','?1','asp_temp_u','humidity_u','##','##2','##3']
            unitsx=['time','seconds_since_1990','pressure_l','Asp_temp_L','humidity_l','WindSpeed_L','winddirection_s_l','swupper','swlower',
                 'lwupper','lwlower','temperatureradsensor','sr_l','sr_u','solar?',
                 'thermistorstring_1','thermistorstring_2','thermistorstring_3','thermistorstring_4','thermistorstring_5','thermistorstring_6','thermistorstring_7','thermistorstring_8',
                 'roll','pitch','heading','mm','gpstime','Latitude','Longitude','Altitude','giodal','geounit?',
                 'battvolt','?1','asp_temp_u','humidity_u','##','##2','##3']
        
            df=pd.read_csv(fn,header=None,names=cols,skiprows=skip)
            df.Latitude=df.Latitude.astype(float)
        
        df[df=="NAN"]=np.nan
        
        if site!='NUK_U':
            for col in cols[2:]:
                df[col] = pd.to_numeric(df[col])
        
        if ((site=='JAR')or(site=='SWC')or(site=='NUK_U')or(site=='QAS_U')):
            sensor_levels=['L']
        else:
            sensor_levels=['L','U']
        
        # if ((site!='NUK_U')or(site!='QAS_U')):
        #     for sensor_level in sensor_levels:
        #         df['kk']=df['Asp_temp_'+sensor_level].astype(float)+273.15
        #         df['tr']=(df['kk']/273.15)**0.5
        #         # tr=tr**0.5
                # df["SR_"+sensor_level]*=df['tr']
        
        # elev=np.nanmean(df.Altitude.astype(float)[:-100])
        
        if ((site!='NUK_U')&(site!='QAS_U')):
            df["date"]=pd.to_datetime(df.time)
            
            df['year'] = df['date'].dt.year
            df['month'] = df['date'].dt.month
            df['day'] = df['date'].dt.day
            df['hour'] = df['date'].dt.hour
                
            df["time"]=pd.to_datetime(df[['year', 'month', 'day', 'hour']])
            df.index = pd.to_datetime(df.time)
        
        # df = df.loc[df['time']<'2021-09-01',:] 

        ##%% lat lon slow
        df['Lat_decimal']=np.nan
        df['Lon_decimal']=np.nan
        
        df.Latitude[df.Latitude=='nan']=np.nan
        df.Longitude[df.Longitude=='nan']=np.nan
        
        for i in range(len(df)):
            if np.isfinite(df.Latitude[i]):
                lat_min=(df.Latitude[i]/100-int(df.Latitude[i]/100))*100
                lon_min=(df.Longitude[i]/100-int(df.Longitude[i]/100))*100
                df['Lat_decimal'][i]=int(df.Latitude[i]/100)+lat_min/60
                df['Lon_decimal'][i]=int(df.Longitude[i]/100)+lon_min/60
        
        #----- compute time dependence of position
        
        i_year=2020
        f_year=2022
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
                v=((df.year==yy+i_year)&(df.month==mm+1)&(np.isfinite(df['Lat_decimal'])))
                means[0,yy,mm]=np.nanmean(df['Lat_decimal'][v])
                counts[0,yy,mm]=np.sum(v)
                v=((df.year==yy+i_year)&(df.month==mm+1)&(np.isfinite(df['Lon_decimal'])))
                means[1,yy,mm]=np.nanmean(df['Lon_decimal'][v])
                counts[1,yy,mm]=np.sum(v)
                v=((df.year==yy+i_year)&(df.month==mm+1)&(np.isfinite(df['Altitude'])))
                means[2,yy,mm]=np.nanmean(df['Altitude'][v])
                counts[2,yy,mm]=np.sum(v)
                means[3,yy,mm]=np.nanstd(df['Altitude'][v])
                print(yy+i_year,mm,means[0,yy,mm],means[1,yy,mm],means[2,yy,mm])
                

                years[cc]=yy+i_year
                months[cc]=mm+1
                dec_year[cc]=yy+i_year+(mm/12)
                lat_count[cc]=counts[0,yy,mm]
                lon_count[cc]=counts[1,yy,mm]
                elev_count[cc]=counts[2,yy,mm]
                thresh=0.8
                if site=='JAR':thresh=0.15
                if lat_count[cc]>744*thresh:
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

                if lon_count[cc]>744*thresh:lon_mean[cc]=means[1,yy,mm]
                if elev_count[cc]>744*thresh:elev_mean[cc]=means[2,yy,mm]
                if elev_count[cc]>744*thresh:elev_std[cc]=means[3,yy,mm]
 
        
                cc+=1
        
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
        df1d['elev'] = df1d['elev'].apply(lambda x: '%.1f' % x)
        df1d["elev std"]=pd.Series(elev_std)
        df1d['elev std'] = df1d['elev std'].apply(lambda x: '%.1f' % x)
        df1d["elev count"]=pd.Series(elev_count)
        ofile='./output/'+site+'_positions_monthly'
        df1d.to_csv(ofile+'.csv',index=None)
        print(df1d)
        # df1d.to_excel(ofile+'.xlsx')
        
        # v=df.Latitude[np.isfinite(df.Latitude)]
        # lat=np.nanmean(df.Lat_decimal[:-24])
        # lon=-np.nanmean(df.Lon_decimal[:-24])

        # elevx=np.nanmean(df.Altitude[:-24])
        
        # print(lat,lon,elev)

        

#%%

# graphics layout option
ly='p' # p for .png, x for console only
wo=1 # write out stats
write_fig=1 # write out figures

sites=['SWC','SDM','NSE','CP1','NAU','CEN','NEM','JAR']
# sites=['SWC']

n_sites=len(sites)

stats=np.zeros((7,n_sites))
date0=['']*n_sites
date1=['']*n_sites
years=np.zeros(n_sites)

for site_index,site in enumerate(sites):
    if site_index!=10:
#     # if site=='SWC':
    # if site=='CP1':
    # if site=='SDM':

        df=pd.read_csv('./output/'+site+'_positions_monthly.csv')
        for col in df.columns:
            df[col] = pd.to_numeric(df[col])
        df['day']=15
        df['date']=pd.to_datetime(df[["year", "month","day"]])
        
        df['elev'][df['elev']==0]=np.nan

        v=np.where(~np.isnan(df['lat']))
        v=v[0]
        n=len(v)
        print(v[0],df['date'][v[0]])
        x=[df['date'][v[0]],df['date'][v[n-1]]]

        lat=df.lat
        lon=df.lon
        elev=df.elev
        st=site_index

        plt.close()
        fig, ax = plt.subplots(3,1,figsize=(10,15))


        cc=0
        stnam=site
        if ly!='n':
            ax[cc].plot(df['date'],df['lat'],'.',color='k')
            ax[0].set_title(stnam+' latitude')
            ax[cc].set_xticklabels([])

            cc+=1
            ax[cc].plot(df['date'],df['lon'],'.',color='k')
            ax[cc].set_title(stnam+' longitude')
            ax[cc].set_xticklabels([])

            cc+=1
            ax[cc].plot(df['date'],df['elev'],'.',color='k')
            ax[cc].set_ylabel('meters a.s.l')
        

        cc=0        
        stats[0,st]=lat[v[0]] 
        ax[cc].plot(df['date'][v[0]],stats[0,st],'s',color='g',alpha=0.5,label='first valid datum')        
        stats[1,st]=lat[v[n-1]]
        ax[cc].plot(df['date'][v[n-1]],stats[1,st],'s',color='r',alpha=0.5,label='last valid datum')
        y=[stats[0,st],stats[1,st]]
        ax[cc].plot(x, y, '--',c='b',label='approximation')
        ax[cc].legend()

        cc+=1

        stats[2,st]=lon[v[0]] 
        ax[cc].plot(df['date'][v[0]],stats[2,st],'s',color='g',alpha=0.5)        
        stats[3,st]=lon[v[n-1]]
        ax[cc].plot(df['date'][v[n-1]],stats[3,st],'s',color='r',alpha=0.5)

        y=[stats[2,st],stats[3,st]]
        ax[cc].plot(x, y, '--',c='b')

        cc+=1

        stats[4,st]=elev[v[0]] 
        ax[cc].plot(df['date'][v[0]],stats[4,st],'s',color='g',alpha=0.5)        
        stats[5,st]=elev[v[n-1]]
        ax[cc].plot(df['date'][v[n-1]],stats[5,st],'s',color='r',alpha=0.5)
        
        plt.setp(ax[cc].xaxis.get_majorticklabels(), rotation=90,ha='center' )
        # ax[cc].xaxis.set_major_formatter(mdates.DateFormatter('%Y %b %d'))

        de=stats[5,st]-stats[4,st]

        ax[2].set_title(stnam+' elevation, change = '+str('%.0f'%de)+' m')

        y=[stats[4,st],stats[5,st]]
        ax[cc].plot(x, y, '--',c='b')
        
        coords_1 = (stats[0,st],stats[2,st])
        coords_2 = (stats[1,st],stats[3,st])
        stats[6,st]=geopy.distance.distance(coords_1, coords_2).m
        # print(df.name[k],df.name[k],str("%.0f"%(dist)))
        
        date0[st]=df['date'][v[0]]
        date1[st]=df['date'][v[n-1]]
        
        duration = df['date'][v[n-1]] - df['date'][v[0]]                         # For build-in functions
        duration_in_s = duration.total_seconds()
        years[st] = divmod(duration_in_s, 86400)[0]/365.25
        print(years[st])


        if ((ly=='p')&(write_fig)): plt.savefig('./figs/'+stnam+'.png', dpi=72,bbox_inches = 'tight')
        if ly=='x':plt.show()

#%% output multi-site table

df2 = pd.DataFrame(columns=['site name','first valid date','latest valid date','delta time','first valid latitude, °N','latest valid latitude, °N','first valid longitude, °W','latest valid longitude, °W','displacement, m','displacement rate, m/y','first valid elevation, m','latest valid elevation, m','elevation change, m'])
df2['site name']=pd.Series(sites)
df2['first valid date']=pd.Series(date0[:])
df2['latest valid date']=pd.Series(date1[:])
df2['first valid latitude, °N']=pd.Series(stats[0,:])
df2['latest valid latitude, °N']=pd.Series(stats[1,:])
df2['first valid longitude, °W']=pd.Series(stats[2,:])
df2['latest valid longitude, °W']=pd.Series(stats[3,:])
df2['first valid elevation, m']=pd.Series(stats[4,:])
df2['latest valid elevation, m']=pd.Series(stats[5,:])
df2["displacement, m"]=pd.Series(stats[6,:])
df2['elevation change, m']=df2['latest valid elevation, m']-df2['first valid elevation, m']
df2["delta time"]=pd.Series(years[:])

df2["displacement rate, m/y"]=df2["displacement, m"]/df2["delta time"]
df2["displacement rate, m/y"][df2['site name']=='KAN_B']=np.nan
df2['displacement, m'] = df2['displacement, m'].map(lambda x: '%.0f' % x)
df2['elevation change, m'] = df2['elevation change, m'].map(lambda x: '%.0f' % x)
df2['delta time'] = df2['delta time'].map(lambda x: '%.1f' % x)
df2['first valid elevation, m'] = df2['first valid elevation, m'].map(lambda x: '%.0f' % x)
df2['latest valid elevation, m'] = df2['latest valid elevation, m'].map(lambda x: '%.0f' % x)
df2['first valid latitude, °N'] = df2['first valid latitude, °N'].map(lambda x: '%.6f' % x)
df2['latest valid latitude, °N'] = df2['latest valid latitude, °N'].map(lambda x: '%.6f' % x)
df2['first valid longitude, °W'] = df2['first valid longitude, °W'].map(lambda x: '%.6f' % x)
df2['latest valid longitude, °W'] = df2['latest valid longitude, °W'].map(lambda x: '%.6f' % x)
df2['displacement rate, m/y'] = df2['displacement rate, m/y'].map(lambda x: '%.1f' % x)

if wo:df2.to_csv('./output/GCN_positions_distance_stats.csv',sep=';')
if wo:df2.to_excel('./output/GCN_positions_distance_stats.xlsx')

print("average displacement rate",np.nanmean(df2["displacement rate, m/y"].astype(float)))
print("average displacement rate",np.nanstd(df2["displacement rate, m/y"].astype(float)))
print("elevation change, m",np.nanmean(df2["elevation change, m"].astype(float)))
print("elevation change, m",np.nanstd(df2["elevation change, m"].astype(float)))

