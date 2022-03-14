#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: jeb@geus.dk


"""

sites=['SWC','CP1','SDM','NSE','NUK_U']

site='CP1'
site='SDM'
# site='NSE'
# site='SWC'
# site='NUK_U'

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

# -------------------------------- chdir
if os.getlogin() == 'adrien':
    base_path = '/home/adrien/EO-IO/rain_optics_SICE_AWS'
elif os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/rain_optics_SICE_AWS'

os.chdir(base_path)

sys.path.append(base_path)
import WetBulb as wb


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

sites=['SWC','SDM','NSE','CP1','NAU','CEN']

for site_index,site in enumerate(sites):
    # if site_index!=10:
#     # if site=='SWC':
    # if site=='CP1':
    if site=='SDM':
    # if site=='CEN':

        if site=='SWC':
            a_or_b='(a)'
            site2='Swiss Camp'
            IMEI='300534060524310' ; lat=69+33.21479/60 ;lon=-(49+22.297377/60); elev=1148
            swd_coef=12.00 ; swu_coef=18.00
            t0_pre_rain=datetime(2021, 8, 9) ; t1_pre_rain=datetime(2021, 8, 11)
            t0_post_rain=datetime(2021, 8, 15) ; t1_post_rain=datetime(2021, 8, 19)
            t0=datetime(2021, 8, 10) ; t1=datetime(2021, 8, 24)
        
        if site=='NUK_U':
            a_or_b='(a)'
            site2='NUK_U'
            xtit='UTC time, year 2021'
            t0=datetime(2021, 8, 9,23) ; t1=datetime(2021, 8, 21)
            t0=datetime(2021, 8, 7) ; t1=datetime(2021, 8, 22)
            t0_pre_rain=datetime(2021, 8, 9) ; t1_pre_rain=datetime(2021, 8, 11)
            t0_post_rain=datetime(2021, 8, 15) ; t1_post_rain=datetime(2021, 8, 19)        
            swd_coef=1 ; swu_coef=1
            
        if site=='SDM':
            a_or_b='(b)'
            IMEI='300534060529220' ; site2='South Dome' ; lat=63.14889	;lon=-44.81667; elev=2895
            t0_pre_rain=datetime(2021, 8, 11) ; t1_pre_rain=datetime(2021, 8, 12)
            t0_post_rain=datetime(2021, 8, 15) ; t1_post_rain=datetime(2021, 8, 19)
            t0=datetime(2021, 8, 11,6) ; t1=datetime(2021, 8, 20,18)
            t0=datetime(2021, 8, 6,12) ; t1=datetime(2021, 8, 30,18)
            t0=datetime(2021, 6, 21,20) ; t1=datetime(2021, 8, 30,18)
            swd_coef=14.75 ; swu_coef=14.44
        
        if site=='NSE':
            a_or_b='(b)'
            IMEI='300534062416350' ; site2='NASA-SE' ; lat=66.2866629	;lon=-42.2961828; elev=2386
            swd_coef=12.66 ; swu_coef=13.9
            t0_pre_rain=datetime(2021, 8, 12,12) ; t1_pre_rain=datetime(2021, 8, 13,12)
            t0_post_rain=datetime(2021, 8, 15,12) ; t1_post_rain=datetime(2021, 8, 16)
            t0=datetime(2021, 8, 13,6) ; t1=datetime(2021, 8, 16,11)
            swd_coef=12.66 ; swu_coef=13.90
        
        if site=='CP1':
            a_or_b='(b)'
            IMEI='300534062418810' ; site2='Crawford Pt. (CP1)'; lat=69.8819; lon=-46.97358; elev=1958
            t0_pre_rain=datetime(2021, 8, 11) ; t1_pre_rain=datetime(2021, 8, 12)
            t0_post_rain=datetime(2021, 8, 15) ; t1_post_rain=datetime(2021, 8, 17)
            t0=datetime(2021, 8, 13) ; t1=datetime(2021, 8, 16)
            t0=datetime(2021, 8, 10) ; t1=datetime(2021, 8, 24)
            swd_coef=15.72 ; swu_coef=12.18
                
        
        if site=='CEN':
            IMEI='300534062419950'; site2='Camp Century' #; lat=69.8819; lon=-46.97358; elev=1958
            t0=datetime(2021, 8, 12,21) ; t1=datetime(2021, 8, 16,18)
            t0=datetime(2021, 8, 12,23) ; t1=datetime(2021, 8, 20)
            t0_pre_rain=datetime(2021, 8, 13) ; t1_pre_rain=datetime(2021, 8, 14)
            t0_post_rain=datetime(2021, 8, 15) ; t1_post_rain=datetime(2021, 8, 17)
        
        if site=='NAU':
            IMEI='300534062412950'; site2='NASA-U' #; lat=69.8819; lon=-46.97358; elev=1958
            t0=datetime(2021, 8, 12,21) ; t1=datetime(2021, 8, 16,18)
            t0_pre_rain=datetime(2021, 8, 18) ; t1_pre_rain=datetime(2021, 8, 19)
            t0_post_rain=datetime(2021, 8, 18) ; t1_post_rain=datetime(2021, 8, 20)
            t0=datetime(2021, 8, 18) ; t1=datetime(2021, 11, 10)
        
        # if site=='NEM':
        #     IMEI='300534062417910' ; site2='NEEM' ; lat=77.5022		;lon=-50.8744; elev=2454
        #     t0=datetime(2021, 8, 13,14) ; t1=datetime(2021, 8, 18,12)
        #     t0_pre_rain=datetime(2021, 8, 13) ; t1_pre_rain=datetime(2021, 8, 14)
        #     t0_post_rain=datetime(2021, 8, 15) ; t1_post_rain=datetime(2021, 8, 17)
        #     t0=datetime(2021, 8, 14) ; t1=datetime(2021, 8, 20)
        
        # if site=='QAS_U':
        #     site2='QAS_U'
        #     xtit='UTC time, year 2021'
        
        # -------------- copy of raw file names to an easier to understand format
        if ((site!='NUK_U')&(site!='QAS_U')):
            fn='/Users/jason/Dropbox/AWS/GCNv2_xmit/AWS_'+IMEI+'.txt' 
            fn='/Users/jason/Dropbox/AWS/xmit_latest/AWS_'+IMEI+'.txt' 
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
        
            df.rename({'AirPressure(hPa)': 'Pressure_L'}, axis=1, inplace=True)
            df.rename({'RelativeHumidity(%)': 'Humidity_L'}, axis=1, inplace=True)
        
            df.rename({'HeightStakes(m)': 'SR_L'}, axis=1, inplace=True)
            df.rename({'HeightSensorBoom(m)': 'SR_U'}, axis=1, inplace=True)
            df.rename({'AirTemperature(C)': 'Asp_temp_L'}, axis=1, inplace=True)
            df.rename({'Ac_Rain(mm)': 'Rain_amount_L'}, axis=1, inplace=True)
        
            df.rename({'LongwaveRadiationDown(W/m2)': 'LWD'}, axis=1, inplace=True)
            df.rename({'LongwaveRadiationUp(W/m2)': 'LWU'}, axis=1, inplace=True)
            df.rename({'ShortwaveRadiationDown_Cor(W/m2)': 'SW Downward'}, axis=1, inplace=True)
            df.rename({'ShortwaveRadiationUp_Cor(W/m2)': 'SW Upward'}, axis=1, inplace=True)
            df.rename({'ElevationGPS(m)': 'Altitude'}, axis=1, inplace=True)
            df.rename({'LongitudeGPS(degW)': 'Longitude'}, axis=1, inplace=True)
            df.rename({'LatitudeGPS(degN)': 'Latitude'}, axis=1, inplace=True)
            df.rename({'WindSpeed(m/s)': 'WindSpeed_L'}, axis=1, inplace=True)
        
            df.rename({'Year': 'year'}, axis=1, inplace=True)
            df.rename({'MonthOfYear': 'month'}, axis=1, inplace=True)
            df.rename({'DayOfMonth': 'day'}, axis=1, inplace=True)
            df.rename({'HourOfDay(UTC)': 'hour'}, axis=1, inplace=True)
        
            df["time"]=pd.to_datetime(df[['year', 'month', 'day', 'hour']])
            df.index = pd.to_datetime(df.time)
        
            # print(df.columns)
            # print(df)
            
        if site=='SWC':
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
        
        if ((site=='SWC')or(site=='NUK_U')or(site=='QAS_U')):
            sensor_levels=['L']
        else:
            sensor_levels=['L','U']
        
        if ((site!='NUK_U')or(site!='QAS_U')):
            for sensor_level in sensor_levels:
                df['kk']=df['Asp_temp_'+sensor_level].astype(float)+273.15
                df['tr']=(df['kk']/273.15)**0.5
                # tr=tr**0.5
                df["SR_"+sensor_level]*=df['tr']
        
        elev=np.nanmedian(df.Altitude.astype(float)[:-100])
        
        if ((site!='NUK_U')&(site!='QAS_U')):
            df["date"]=pd.to_datetime(df.time)
            
            df['year'] = df['date'].dt.year
            df['month'] = df['date'].dt.month
            df['day'] = df['date'].dt.day
            df['hour'] = df['date'].dt.hour
                
            df["time"]=pd.to_datetime(df[['year', 'month', 'day', 'hour']])
            df.index = pd.to_datetime(df.time)
        
        # df = df.loc[df['time']>'2021-08-01',:] 
        # df = df.loc[df['time']<'2021-09-01',:] 
        
        ##%% albedo
        N=len(df)
        # ---------------------------------------------- compute albedo
        df['ALB']=np.nan
        plusminus=11
        df['SW Downward']*=10
        df['SW Upward']*=10
        df['SRnet']=df['SW Downward']-df['SW Upward']
        df['SW Downward'][df['SRnet']<0]=np.nan
        df['SW Upward'][df['SRnet']<0]=np.nan
        # plt.plot(df['SW Downward']-df['SW Upward'])
        # plt.plot(df['SW Downward']-df['SW Upward'])
        
        for i in range(0+plusminus,N-plusminus):
            df['ALB'][i]=np.nansum(df['SW Upward'][i-plusminus:i+plusminus])/np.nansum(df['SW Downward'][i-plusminus:i+plusminus])
        
        # df['ALB'][df['ALB']>0.95]=np.nan
        # df['ALB'][df['ALB']<0.3]=np.nan
        
        ##%% lat lon slow
        df['Lat_decimal']=np.nan
        df['Lon_decimal']=np.nan
        
        df.Latitude[df.Latitude=='nan']=np.nan
        df.Longitude[df.Longitude=='nan']=np.nan
        
        
        for i in range(N):
            if np.isfinite(df.Latitude[i]):
                lat_min=(df.Latitude[i]/100-int(df.Latitude[i]/100))*100
                lon_min=(df.Longitude[i]/100-int(df.Longitude[i]/100))*100
                df['Lat_decimal'][i]=int(df.Latitude[i]/100)+lat_min/60
                df['Lon_decimal'][i]=int(df.Longitude[i]/100)+lon_min/60
        
        v=df.Latitude[np.isfinite(df.Latitude)]
        lat=np.nanmedian(df.Lat_decimal[:-24])
        lon=-np.nanmedian(df.Lon_decimal[:-24])

        elevx=np.nanmedian(df.Altitude[:-24])
        
        print(lat,lon,elev)
        
#         #%%
#         opath='/Users/jason/Dropbox/AWS/GCNet_positions/output/'
#         out=open(opath+site+'_latest_position.csv','w+')
#         out.write('site,nickname,latest date,latitude decimal,longitude decimal,latitude degrees,latitude decimal minutes,longitude degrees,longitude decimal minutes,elevation\n')
#         out.write(site2+','+\
#             site+','+\
#               df.time[-1].strftime('%Y-%m-%d')+','+\
#               "{:.4f}".format(lat)+','+\
#               "{:.4f}".format(lon)+','+\
#               "{:.0f}".format(int(lat))+','+\
#               "{:.4f}".format(lat_min)+','+\
#               "{:.0f}".format(int(lon))+','+\
#               "{:.4f}".format(lon_min)+','+\
#               "{:.0f}".format(elevx)
#               )
#         out.close()
#         #%%

#         kml = simplekml.Kml(open=1)
    
#         pnt = kml.newpoint(name=site)
#         pnt.coords=[(lon,lat)]
#         pnt.style.labelstyle.color = simplekml.Color.red  # Make the text red
#         pnt.style.labelstyle.scale = 1.5  # text scaling multiplier
#         pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
#         pnt.altitudemode = simplekml.AltitudeMode.relativetoground

#         kml_ofile=opath+site+".kml"
#         kml.save(kml_ofile)

# #%% read instrument height data
# from glob import glob

# files = sorted(glob(opath+"*.csv"), reverse=False)

# z_df=pd.read_csv(files[0])
# for i,file in enumerate(files):
#     print(i,file)
#     if i>0:
#         z_df = pd.concat([z_df,pd.read_csv(files[i])])

# print(z_df)

# z_df.to_excel(opath+'GCN_latest_positions.xlsx',index=None)
