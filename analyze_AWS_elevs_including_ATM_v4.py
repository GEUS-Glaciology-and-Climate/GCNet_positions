#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 09:51:29 2024

@author: jason
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 07:32:22 2023

@author: jason, jeb@geus.dk

compare ATM and other sources of AWS elevations

see ATM folder in this repository

for the GEUS GPS, this scripts reads output from which currently has some simple filtering of SDL and CP1 outliers
    ./GCN_positions_timeseries_v_thredds.py

"""

from glob import glob
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
from pathlib import Path
import calendar
import simplekml

def save_kml(lat,lon,namex,opath,ofile):
    kml = simplekml.Kml(open=1)
    pnt = kml.newpoint(name=namex)
    pnt.coords=[(lon,lat)]
    pnt.style.labelstyle.color = simplekml.Color.blue  # Make the text red
    pnt.style.labelstyle.scale = 1.5  # text scaling multiplier
    pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
    pnt.altitudemode = simplekml.AltitudeMode.relativetoground
    
    kml_ofile=opath+ofile+".kml"
    kml.save(kml_ofile)


# ## change to your system's login name to change dir for local work
if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/'

os.chdir(base_path)

# read Greenland AWS locations
meta = pd.read_csv('/Users/jason/Dropbox/AWS/GCNET/GC-Net-level-1-data-processing/L1/GC-Net_location.csv')
meta = meta.rename({'Name': 'name'}, axis=1)
meta = meta.rename({'Latitude (°N)': 'lat'}, axis=1)
meta = meta.rename({'Longitude (°E)': 'lon'}, axis=1)
meta = meta.rename({'Elevation (wgs84 m)': 'elev'}, axis=1)

meta["nickname"]=''

sites=['Swiss Camp', 'Crawford Point 1', 'CP2', 'JAR1', 'JAR2', 'JAR3',
'NASA-U', 'GITS', 'Humboldt', 'Summit', 'Tunu-N', 'DYE-2',
'Saddle', 'South Dome', 'NASA-E', 'NASA-SE', 'NGRIP', 'NEEM',
'EastGRIP', 'KAR', 'KULU', 'Aurora', 'Petermann Glacier',
'Petermann ELA']

nicknames=['SWC', 'CP1', 'CP2', 'JAR', 'JR2', 'JR3',
'NAU', 'GIT', 'HUM', 'SUM', 'TUN', 'DY2',
'SDL', 'SDM', 'NAE', 'NSE', 'NGRP', 'NEM',
'EGP', 'KAR', 'KUL', 'AUR', 'PTG',
'PET']

for ss,site in enumerate(sites):
    meta["nickname"][meta["name"] == sites[ss]]=nicknames[ss]

names=['SMS-PET','SMS1','SMS2','SMS3','SMS4','SMS5','LAR1','LAR2','LAR3','Swiss Camp 10m','Roof_GEUS','LYN_T','LYN_L','CEN1','CEN2','KPC_Lv3','KPC_Uv3','THU_L2','WEG_B','ZAK_Uv3']
for name in names:
    meta.drop(meta[meta.name==name].index, inplace=True) # drop original time column

print(meta.name)

meta['elev_linear_slope']=np.nan
meta['elev_linear_intercept']=np.nan
meta['elev_fit_t0']=np.nan
meta['elev_fit_t1']=np.nan
meta['elev_change_linear']=np.nan
meta['elev_change_n_years']=np.nan
meta['elev_mean_from_altimetry']=np.nan
meta['k']=np.nan
meta['N_altimetry_measurements']=np.nan

print(meta.columns)
#%%



# obtain geoidal heights using data from  https://www.agisoft.com/downloads/geoids/
# EGM = earth gravitational model
import rasterio
geoids=['egm2008_25','egm96_15']
geoids=['egm96_15']
# geoids=['egm2008_25']
# geoids=['egm2008_25']

for geoid in geoids:
    dat = rasterio.open(f"/Users/jason/0_dat/geoid/us_nga_{geoid}.tiff")
    # read all the data from the first band
    z = dat.read()[0]
    
    #%%
    # check the crs of the data
    # dat.crs
    # >>> CRS.from_epsg(4326)
    
    # check the bounding-box of the data
    # dat.bounds
    # >>> Out[49]: BoundingBox(left=-120.0, bottom=45.0, right=-117.0, top=48.0)
    
    # since the raster is in regular lon/lat grid (4326) we can use 
    # `dat.index()` to identify the index of a given lon/lat pair
    # (e.g. it expects coordinates in the native crs of the data)
    
    def getval(lon, lat):
        idx = dat.index(lon, lat, precision=1E-6)    
        # return dat.xy(*idx), z[idx]
        return z[idx]
    
    
    N=len(meta)
    x=np.zeros(N)
    
    for i in range(N):
       # print( getval(meta.lon.values[i], meta.lat.values[i]))
       x[i]=getval(meta.lon.values[i], meta.lat.values[i])
       # print(i,meta.lon.values[i], meta.lat.values[i],x)
    
    meta[geoid]=x
    
    #%%
    df = pd.read_excel('./meta/GC-Net historical positions.xlsx')
    os.system('open '+'./meta/GC-Net historical positions.xlsx')
    df.drop(df[df.date == 'doc_2000'].index, inplace=True)
    df.drop(df[df.date == 'ref'].index, inplace=True)
    df.drop(df[df.date == 'WMO-DMI_2012'].index, inplace=True)
    df.drop(df[df.date == 35596].index, inplace=True)
    df.drop(df[df.elev == '-'].index, inplace=True)
    df.drop(df[df.elev == np.nan].index, inplace=True)
    
    df.loc[ df["name_long"] == "Crawford Pt. 1","name_long"] = "Crawford Point 1"
    
    df.date=pd.to_datetime(df.date)
    df['year']=df.date.dt.year
    df['doy']=df.date.dt.dayofyear
    df['jy']=df['year']+df['doy']/366
    # print(np.array(df.date))
    # print(len(df))
    # df.columns
    # v=np.where(((df.name_long=='JAR1')&(df.elev>940)))
    # df.elev[v[0]]=np.nan
    vandecrux_position_compilation=df.copy()
    #%%
    
    fn='/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/data/GNSS data/2023/Data_overview_2023.xlsx'
    GNSS_2023 = pd.read_excel(fn,sheet_name='Calculated_positions',skiprows=1)
    # GNSS_2023_meta = pd.read_excel(fn,sheet_name='Time_and_meta')
    GNSS_2023 = GNSS_2023.rename({'Station name': 'name'}, axis=1)
    GNSS_2023 = GNSS_2023.rename({'Ortho H': 'elev'}, axis=1)
    # GNSS_2023_meta['date']=pd.to_datetime(GNSS_2023_meta['Date'])
    # print(GNSS_2023_meta.columns)
    GNSS_2023['date']=pd.to_datetime('20230615',format='%Y%m%d')
    GNSS_2023['year']=GNSS_2023.date.dt.year
    GNSS_2023['doy']=GNSS_2023.date.dt.dayofyear
    GNSS_2023['jy']=GNSS_2023['year']+GNSS_2023['doy']/366
    

    print(GNSS_2023.columns)
    #%%
    n_AWS=len(meta)
    
    # meta=meta.sort_values(by='name', ascending=True)
    # print(meta.columns)
    
    
    # iyear=1995 ; fyear=1998
    # n_years=fyear-iyear+1
    # years=np.arange(iyear,fyear+1).astype(str)
    
    years=['1995',
    '1996',
    '1997',
    '1998',
    '1999',
    '2001',
    '2002',
    '2005',
    '2006',
    '2007',
    '2008',
    '2010',
           '2011',
           '2012',
           '2013',
           '2014',
           '2015',
           '2016',
           '2017',
           '2018',
           '2019',
    ]
    
    n_years=len(years)
    
    # print(n_years)
    
    #%%
    
    th=1 
    font_size=10
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
    ms=12
    
    sites=['Swiss Camp', 'Crawford Point 1', 'CP2', 'JAR1', 'JAR2', 'JAR3',
    'NASA-U', 'GITS', 'Humboldt', 'Summit', 'Tunu-N', 'DYE-2',
    'Saddle', 'South Dome', 'NASA-E', 'NASA-SE', 'NGRIP', 'NEEM',
    'EastGRIP', 'KAR', 'KULU', 'Aurora', 'Petermann Glacier',
    'Petermann ELA']

    # site names compatible with interpolated lat lons
    sites2=['Swiss Camp', 'Crawford Pt. 1', 'CP2', 'JAR1', 'JAR2', 'JAR3',
    'NASA-U', 'GITS', 'Humboldt', 'Summit', 'Tunu-N', 'DYE-2',
    'Saddle', 'South Dome', 'NASA-E', 'NASA-SE', 'NGRIP', 'NEEM',
    'East Grip', 'KAR', 'KULU', 'Aurora', 'Petermann Glacier',
    'Petermann ELA']    
    
    constrain_elev=0
    do_GEUS_GPS=1
    show_approximation_elev_history=1
    show_khan=0 # switches to 1 if SWC or JAR
    show_stober=0 # switches to 1 if SWC or JAR
    show_ATM=1
    

    choice_site=0 ; elev0=1116 ; elev1=1157 # SWC
    choice_site=1 # CP1
    choice_site=2 # CP2
    # choice_site=3 # JAR1
    # choice_site=4 # JAR2
    # choice_site=5 # JAR3
    choice_site=6 # NAU
    # choice_site=7 # GIT
    choice_site=8 # HUM
    # choice_site=9 # SUM
    choice_site=10 # TUN
    # choice_site=11 # DY2
    # choice_site=12 # SDL
    # choice_site=13 # SDM
    # choice_site=14 # NAE
    # choice_site=15 # NSE
    # choice_site=16 # NGP
    # choice_site=17 # NEM
    # choice_site=18 # EGP
    # choice_site=19 # KAR
    # choice_site=20 # KUL
    # choice_site=21 # AUR
    # choice_site=22 # PTG
    # choice_site=23 # PTE
    
    elevs=np.zeros(n_years)
    time_ATM_decimal_year=np.zeros(n_years)
    time_ATM_decimal_year_v2=np.zeros(n_years)
    dist=np.zeros(n_years)

    def elev_v_time_interpolation(fn,sites2,k):
        print(f'elev_v_time_interpolation ',sites2[k])
        pos_v_time=pd.read_csv(fn)
        # implement elevation approximation
        n_time_epochs=len(time_elev_approximation)
        pos_v_time['date']=pd.to_datetime(pos_v_time['date'])
        pos_v_time['year']=pos_v_time['date'].dt.year
        pos_v_time['doy']=pos_v_time['date'].dt.day_of_year
        pos_v_time['n_days']=365
        for m in range(len(pos_v_time)):
            if calendar.isleap(pos_v_time.year[m]):
                pos_v_time['n_days']=366
        pos_v_time['jdy']=pos_v_time['year']+pos_v_time['doy']/pos_v_time['n_days']
        
        pos_v_time['elev']=np.nan
        for l in range(n_time_epochs-1):
            # print(l)
            dz=elev_approximation[l+1]-elev_approximation[l]
            dx=time_elev_approximation[l+1]-time_elev_approximation[l]
            m=dz/dx
            b=elev_approximation[l]-m*time_elev_approximation[l]
            v=((pos_v_time['jdy']>=time_elev_approximation[l])&(pos_v_time['jdy']<time_elev_approximation[l+1]))
            pos_v_time['elev'][v]=m*pos_v_time['jdy'][v]+b
        plt.close()
        plt.plot(pos_v_time['elev'])
        # print(pos_v_time)
        # adjust data precision
        vals=['elev']
        for val in vals:
            pos_v_time[val] = pos_v_time[val].map(lambda x: '%.1f' % x)
        pos_v_time.to_csv(f'./output/{sites2[k]}_position_interpolated_with_elev.csv',columns=['date','lon','lat','elev'],index=None)
        return None
        
    # n_AWS=1
    # for k in range(n_AWS):
    for k in [choice_site]:
        
        min_tolerated_dist=1

        if nicknames[k]=='SWC':
            show_stober=1
            show_khan=1
            
        if nicknames[k]=='JAR':
            show_stober=1
            show_khan=1
        
        GNSS_2023_elev_pick=np.nan
        v=nicknames[k]==GNSS_2023.name
        if np.sum(v)>0:
            GNSS_2023_decimal_year=GNSS_2023['jy'][v]
            GNSS_2023_elev_pick=GNSS_2023.elev.values[v][0]

        if nicknames[k]=='SWC':
            time_elev_approximation=[1990.0,1990.5,2005,2015,2020,2023.5]
            elev_approximation=[1155.5,1155.5,1135,1122,1119.2,1119.2]
        if nicknames[k]=='CP1':
            time_elev_approximation=[1995,2023.5]
            elev_approximation=[1960,GNSS_2023_elev_pick]
        if nicknames[k]=='JAR':
            time_elev_approximation=[1996.5,2010,2015,2023.5]
            elev_approximation=[932.3,919.5,909,906]
        if nicknames[k]=='JR2':
            time_elev_approximation=[1999,2020]
            elev_approximation=[523,484]
        if nicknames[k]=='JR3':
            time_elev_approximation=[2000,2020]
            elev_approximation=[300,220]
        if nicknames[k]=='NAU':
            time_elev_approximation=[1995,2023.5]
            elev_approximation=[2338,GNSS_2023_elev_pick]
        if nicknames[k]=='GIT':
            time_elev_approximation=[1995,2022]
            elev_approximation=[1872,1871]
        if nicknames[k]=='HUM':
            time_elev_approximation=[1995,2023.5]
            elev_approximation=[1974,GNSS_2023_elev_pick]
            min_tolerated_dist=8
        if nicknames[k]=='SUM':
            time_elev_approximation=[1996,2023]
            elev_approximation=[3205,3205]
        if nicknames[k]=='TUN':
            time_elev_approximation=[1996,2023]
            elev_approximation=[2074,2078]
        if nicknames[k]=='DY2':
            time_elev_approximation=[1996,2023.5]
            elev_approximation=[2115,GNSS_2023_elev_pick]
        if nicknames[k]=='SDL':
            time_elev_approximation=[1997,2023.5]
            elev_approximation=[2456,GNSS_2023_elev_pick]
        if nicknames[k]=='SDM':
            time_elev_approximation=[1997,2023.5]
            elev_approximation=[2878,GNSS_2023_elev_pick]
        if nicknames[k]=='NAE':
            time_elev_approximation=[1997,2023.5]
            elev_approximation=[2622.7,GNSS_2023_elev_pick]
        if nicknames[k]=='NSE':
            time_elev_approximation=[1998,2023.5]
            elev_approximation=[2372,GNSS_2023_elev_pick]
            min_tolerated_dist=4
        if nicknames[k]=='NGRP':
            time_elev_approximation=[1997,2021]
            elev_approximation=[2919.3,2919.3]
            min_tolerated_dist=4
        if nicknames[k]=='NEM':
            time_elev_approximation=[2006,2023.5]
            elev_approximation=[2451,2450.2]
            min_tolerated_dist=4
        if nicknames[k]=='EGP':
            time_elev_approximation=[2014,2023.5]
            elev_approximation=[2663.5,2663.2]
            min_tolerated_dist=4
        if nicknames[k]=='KAR':
            time_elev_approximation=[1999,2002]
            elev_approximation=[2552.5,2552.5]
            min_tolerated_dist=4
        if nicknames[k]=='KUL':
            time_elev_approximation=[1999,2001]
            elev_approximation=[900,900]
            min_tolerated_dist=2
        if nicknames[k]=='AUR':
            time_elev_approximation=[1999,2001]
            elev_approximation=[1748,1748]
            min_tolerated_dist=5
        if nicknames[k]=='PTG':
            time_elev_approximation=[2002,2007]
            elev_approximation=[33,28]
            min_tolerated_dist=5
        if nicknames[k]=='PET':
            time_elev_approximation=[2002,2021]
            elev_approximation=[949,940]
            min_tolerated_dist=3

        
        # read in interpolated lat lon
        fn=f'./output/{sites2[k]}_position_interpolated.csv'
        elev_v_time_interpolation(fn,sites2,k)

        #%%
        
        if nicknames[k]=='JAR': # JAR
            # v=np.where(((df.name_long=='JAR1')&(df.elev>940)))
            # v=np.where(df.name_long=='JAR1')
            # ix=v[0]
            # for i in ix:
            #     print(pd.to_datetime(df.date.values[i]).strftime('%Y-%m-%d'))
            #     datexx=pd.to_datetime(df.date.values[i]).strftime('%Y-%m-%d')
            #     save_kml(df.lat.values[i],df.lon.values[i],'JAR1_'+datexx,'/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/output/kml/','JAR1_'+datexx)

            geoid_offset=-meta[geoid].values[k]

            fn='./meta/Stober_et_al_2023/SWC+ST2-Alle Pegel_Koordinaten Geodätisch.xlsx'
            # os.system('open '+fn)
            stober=pd.read_excel(fn,skiprows=7,names=['latdeg','latmin','latsec','NS','lat','emt','londeg','lonmin','lonsec','EW','lon','elevation'],index_col=0)
            stober = stober.iloc[65:]
            stober_dates=[]
            stober_elevs=[]
            for i,temp in enumerate(stober.index):
                temp=str(temp)
                stake=temp[0:5]
                # print(i,stake)
                
                # temp.replace('106-ex010915','nan')
                # temp.replace('ex010915','nan')
                # print(i,temp)
                
                if temp!='nan':
                    lat=stober.lat.values[i]
                    lon=-stober.lon.values[i]
                    datex=pd.to_datetime(temp[6:],format=('%d%m%y'))
                    print(stake,datex.strftime('%Y-%m-%d'),lat,lon,stober.elevation.values[i])
                    ST2_offset=35
                    if stake=='ST201':
                        stober_dates.append(datex.strftime('%Y-%m-%d'))
                        # stober_offset=-66
                        stober_elevs.append(stober.elevation.values[i]-1.08+geoid_offset-ST2_offset)
                    # save_kml(lat,lon,stake+' '+datex.strftime('%Y-%m-%d'),'/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/output/kml/Stober_Swiss_Camp/',stake+' '+datex.strftime('%Y-%m-%d'))

            stober2=pd.DataFrame({'date':np.array(stober_dates),
                              'elev':np.array(stober_elevs),
                              })
            stober2['date']=pd.to_datetime(stober2.date)
            stober2['year'] = stober2['date'].dt.year
            stober2['doy'] = stober2['date'].dt.dayofyear
            stober2['n_days']=365
            for m in range(len(stober2)):
                if calendar.isleap(stober2.year[m]):
                    stober2['n_days']=366
            stober2['jdy']=stober2['year']+stober2['doy']/stober2['n_days']

        if nicknames[k]=='SWC': # SWC
            
            fn='./meta/Stober_et_al_2023/SWC+ST2-Alle Pegel_Koordinaten Geodätisch.xlsx'
            # heights are referenced to EUREF, a height system I used since my first measurements in 1991. 
            # For correction to ITRF use ITRF = EUREF  - 1.08 Meter.
            # os.system('open '+fn)
            stober=pd.read_excel(fn,skiprows=7,names=['latdeg','latmin','latsec','NS','lat','emt','londeg','lonmin','lonsec','EW','lon','elevation'],index_col=0)
            stober = stober.iloc[:-32]
            stober_dates=[]
            stober_elev0=[]
            stober_elev1=[]
            stober_elev2=[]
            for i,temp in enumerate(stober.index):
                temp=str(temp)
                pegel=temp[0:3]
                temp.replace('106-ex010915','nan')
                temp.replace('ex010915','nan')
                # print(i)
                if temp!='nan' and temp[4:]!='ex010915':
                    lat=stober.lat.values[i]
                    lon=-stober.lon.values[i]
                    datex=pd.to_datetime(temp[4:],format=('%d%m%y'))
                    # print(pegel,datex.strftime('%Y-%m-%d'),lat,lon,stober.elevation.values[i])
                    pegs=['106','120','121'] 
                    # oss=1.3
                    geoid_offset=-meta[geoid].values[k]
                    geoid_offset=-25
                    for peg in pegs:
                        if pegel=='106':
                            stober_dates.append(datex.strftime('%Y-%m-%d'))
                            # geoid_offset=-25
                            stober_elev0.append(stober.elevation.values[i]-1.08+geoid_offset)
                        if pegel=='120':
                            # geoid_offset=-41/oss
                            stober_elev1.append(stober.elevation.values[i]-1.08+geoid_offset)                    
                        if pegel=='121':
                            # geoid_offset=-44/oss
                            stober_elev2.append(stober.elevation.values[i]-1.08+geoid_offset)
                    # save_kml(lat,lon,pegel+' '+datex.strftime('%Y-%m-%d'),'/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/output/kml/Stober_Swiss_Camp/',pegel+' '+datex.strftime('%Y-%m-%d'))

            stober2=pd.DataFrame({'date':np.array(stober_dates),
                              'elev0':np.array(stober_elev0),
                              'elev1':np.array(stober_elev1),
                              'elev2':np.array(stober_elev2),
                              })
            stober2['date']=pd.to_datetime(stober2.date)
            stober2['year'] = stober2['date'].dt.year
            stober2['doy'] = stober2['date'].dt.dayofyear
            stober2['n_days']=365
            for m in range(len(stober2)):
                if calendar.isleap(stober2.year[m]):
                    stober2['n_days']=366
            stober2['jdy']=stober2['year']+stober2['doy']/stober2['n_days']



        plt.close()
        plt.clf()
        fig, ax = plt.subplots(figsize=(9,7))
        yx=np.zeros(n_years)
        ATM=pd.read_csv(f'./ATM/output/{nicknames[choice_site]}.csv')
    
        ATM.columns
        ATM['date']=pd.to_datetime(ATM.date)
        ATM['doy'] = ATM['date'].dt.dayofyear
        ATM['year'] = ATM['date'].dt.year
        ATM['n_days']=365
        for m in range(len(ATM)):
            if calendar.isleap(ATM.year[m]):
                ATM['n_days']=366
        ATM['jdy']=ATM['year']+ATM['doy']/ATM['n_days']
   
        time_ATM_decimal_year=ATM['jdy'].values
        elevs=ATM.elev_ATM_m-meta[geoid].values[k]
        dist=ATM.distance_km
        
        x1=np.sin(np.radians(ATM.slope_S2N_deg.values))*ATM.distance_km.values*1000
        x2=np.sin(np.radians(ATM.slope_W2E_deg.values))*ATM.distance_km.values*1000
        # print('slope cor',ATM.slope_S2N_deg.values,x1)
        
        # x1=ATM.slope_S2N.values[v][0]*ATM.dist.values[v][0]*1000
        # x2=ATM.slope_W2E.values[v][0]*ATM.dist[v][0]*1000
        x1=ATM.slope_S2N_deg.values*ATM.distance_km.values*1000
        x2=ATM.slope_W2E_deg.values*ATM.distance_km*1000
        yx=ATM.elev_ATM_m.values-meta[geoid].values[k]-x1-x2
        # if dist<min_tolerated_dist:
        #     # print(time_ATM_decimal_year[yy],"%.1f"%elevs[yy],"%.1f"%dist[yy],"%.1f"%(ATM.slope_S2N.values[v][0]*1000),"%.1f"%(ATM.slope_W2E.values[v][0]*1000),
        #     #       "%.1f"%x1,"%.1f"%x2,"%.1f"%(x1+x2))
        #     print(geoid,time_ATM_decimal_year[yy],"%.1f"%elevs[yy],"%.1f"%dist[yy],"%.1f"%(ATM.slope_S2N.values[v]*1000),"%.1f"%(ATM.slope_W2E.values[v]*1000),
        #           "%.1f"%x1,"%.1f"%x2,"%.1f"%(x1+x2))    
        v=np.where(dist<min_tolerated_dist)
        v=v[0]
        
        ATM_dates=time_ATM_decimal_year[v]
        elevs_ATM_within_distance_tolerance=elevs.values[v]
        y2=yx.values[v]
    
        # exclude suspicious values else plot stober
        if nicknames[k]=='SWC':
            # inv=np.where((x<2005)&(y<1125))
            # inv=inv[0]
            # y[inv]=np.nan

            # inv=np.where((x>2015)&(y<1120))
            # inv=inv[0]
            # y[inv]=np.nan
        
            khan=pd.read_csv('/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/Khan/SwissCamp.txt',skiprows=2,delim_whitespace=True,names=['jy','elev'])
            # print(khan)
            # #%%
            if show_stober:
                plt.plot(stober2['jdy'],stober2['elev0'],'s',color='k',label="GPS survey c/o M. Stober, Jan 2024, average of\nstakes 106 with 25m offset: %.1f"%np.mean(stober2['elev0'])+"±%.1f"%np.std(stober2['elev0'])+' m')
            if show_khan:
                os_khan=1135.7
                plt.plot(khan.jy,khan.elev+os_khan,'-',linewidth=th*4,color='orange',label='satellite altimetry c/o S.A. Khan, Jan. 2024\noffset to GEUS GPS data 2022 to 2023',zorder=20)

        # exclude suspicious values else plot stober
        if nicknames[k]=='JAR':
            # inv=np.where((x<2005)&(y<1125))
            # inv=inv[0]
            # y[inv]=np.nan

            # inv=np.where((x>2015)&(y<1120))
            # inv=inv[0]
            # y[inv]=np.nan
            if show_stober:
                plt.plot(stober2['jdy'],stober2['elev'],'s',color='k',zorder=20,
                     label=f"Stober ST201 including {ST2_offset} m offset: %.1f"%np.mean(stober2['elev'])+"±%.1f"%np.std(stober2['elev'])+' m')
            if show_khan:
                khan=pd.read_csv('/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/Khan/JAR.txt',skiprows=2,delim_whitespace=True,names=['jy','elev'])
            # print(khan)
            # #%%
            # plt.plot(stober2['jdy'],(stober2['elev0']+stober2['elev1']+stober2['elev2'])/3,'s',color='k',label="GPS survey c/o M. Stober, Jan 2024, average of\nstakes 106,120,121 with %.1f"%geoid_offset+" m offset: %.1f"%np.mean(stober2['elev0'])+"±%.1f"%np.std(stober2['elev0'])+' m')
            os_khan=927
            plt.plot(khan.jy,khan.elev+os_khan,'-',linewidth=th*4,color='orange',label='satellite altimetry c/o S.A. Khan, Jan. 2024\noffset to GEUS GPS data 2022 to 2023',zorder=20)


        if nicknames[k]=='CP1':
            inv=np.where((ATM_dates>2005)&(elevs_ATM_within_distance_tolerance>1960))
            inv=inv[0]
            elevs_ATM_within_distance_tolerance[inv]=np.nan
        if nicknames[k]=='CP2':
            inv=np.where((ATM_dates>2005)&(elevs_ATM_within_distance_tolerance>1960))
            inv=inv[0]
            elevs_ATM_within_distance_tolerance[inv]=np.nan
        if nicknames[k]=='JR2':
            inv=np.where((ATM_dates>2010)&(elevs_ATM_within_distance_tolerance>530))
            inv=inv[0]
            elevs_ATM_within_distance_tolerance[inv]=np.nan
        if nicknames[k]=='JR3':
            inv=np.where((ATM_dates>2005)&(elevs_ATM_within_distance_tolerance<220))
            inv=inv[0]
            elevs_ATM_within_distance_tolerance[inv]=np.nan
    
        if show_ATM:
            xxx=1000
            lab="ATM within %.1f"%min_tolerated_dist+" km: %.1f"%np.mean(elevs[v])+"±%.1f"%np.std(elevs[v])+' m, using '+geoid.split('_')[0].upper()
            # print('lab',lab)
            # print('elevs',elevs[v])
            plt.plot(ATM_dates,elevs_ATM_within_distance_tolerance,'o', fillstyle='none',markersize=ms, mew=2,label=lab)
            # plt.plot(ATM_dates,y2,'s', fillstyle='none',markersize=ms,label="ATM with slope cor: %.1f"%np.mean(yx[v])+"±%.1f"%np.std(yx[v])+' m')
        
        if ~np.isnan(GNSS_2023_elev_pick):
            lab="Jakobsen et al. GEUS 2023 GNSS: %.1f"%GNSS_2023_elev_pick+'±0.1 m'
            plt.plot(GNSS_2023_decimal_year,GNSS_2023_elev_pick,'o', fillstyle='none',color='c',markersize=ms, mew=2,label=lab,zorder=20)

        # if len(x)>1:
        #     v=np.where(~np.isnan(y))
        #     v=v[0]
        #     x=x[v]
        #     y=y[v]
        #     N_valid=len(y)
        #     b, m = polyfit(x, y, 1)
        #     xx=[x[0],x[-1]]
        #     yy=[xx[0]*m+b,xx[1]*m+b]
        #     dy=yy[1]-yy[0]
        #     ny=xx[1]-xx[0]
        #     plt.plot(xx,yy,c='c',label=f'ATM fit averages {"%.0f"%np.mean(yy)}m')
        #     # print(np.mean(y),np.std(y))
        #     # kx=meta.ID.values[k]-1
        #     # kx=k-1
        #     # kx=np.where()
        #     # kx=kx[0]
        #     # print('hi',meta.name.values[k],k,kx)
        #     meta['elev_linear_slope'][meta.name==sites[k]]=m
        #     meta['elev_linear_intercept'][meta.name==sites[k]]=b
        #     meta['elev_fit_t0'][meta.name==sites[k]]=xx[0]
        #     meta['elev_fit_t1'][meta.name==sites[k]]=xx[1]
        #     meta['elev_change_linear'][meta.name==sites[k]]=dy
        #     meta['elev_change_n_years'][meta.name==sites[k]]=ny
        #     meta['elev_mean_from_altimetry'][meta.name==sites[k]]=np.mean(yy)
        #     meta['N_altimetry_measurements'][meta.name==sites[k]]=N_valid

        # if len(x)==1:
        #     meta['elev_linear_slope'][meta.name==sites[k]]=np.nan
        #     meta['elev_linear_intercept'][meta.name==sites[k]]=np.nan
        #     meta['elev_fit_t0'][meta.name==sites[k]]=np.nan
        #     meta['elev_fit_t1'][meta.name==sites[k]]=np.nan
        #     meta['elev_change_linear'][meta.name==sites[k]]=np.nan
        #     meta['elev_change_n_years'][meta.name==sites[k]]=1
        #     meta['elev_mean_from_altimetry'][meta.name==sites[k]]=y
        #     meta['N_altimetry_measurements'][meta.name==sites[k]]=1

        plt.title(sites[k]+' a.k.a. '+nicknames[k])
        plt.ylabel('AWS barometer elevation above mean sea level, m')
        
        year0=1989 ; year1=2024
        plt.xlim(year0,year1)
        
        v=vandecrux_position_compilation.name_long==sites[k]
        y3=vandecrux_position_compilation.elev[v]
        if ~np.isnan(np.std(y3)):
            # if nicknames[k]=='SWC':
            #     plt.plot(vandecrux_position_compilation.jy[v][0],y3[0],'s', fillstyle='none',markersize=ms/2,c='b',label="Ohmura et al 1991: %.0f"%y3[0]+' m')
            # else:
            plt.plot(vandecrux_position_compilation.jy[v],y3,'s', fillstyle='none',markersize=ms/2,c='b',
                     label="GC-Net historical positions.xlsx: %.0f"%np.mean(y3)+"±%.0f"%np.std(y3)+' m')
        
    
        out=pd.DataFrame({'jdy':ATM_dates,
                          'elev':elevs_ATM_within_distance_tolerance,
                          })
        
        if do_GEUS_GPS:
            suffix=''
            fn = Path(f'./output/Jason/{nicknames[k]}_positions_monthly.csv')
            cat_flag=0
            if fn.is_file():
                # print(fn)
                GEUS_AWS_position=pd.read_csv(fn)
                GEUS_AWS_position.elev[GEUS_AWS_position.elev==0]=np.nan
                GEUS_AWS_position.elev-=1.5 # correction to height of ice surface
                GEUS_AWS_position['day']=15
                
                GEUS_AWS_position['date']=pd.to_datetime(GEUS_AWS_position[['year', 'month', 'day']])
                GEUS_AWS_position['doy'] = GEUS_AWS_position['date'].dt.dayofyear
                GEUS_AWS_position['n_days']=365
                for m in range(len(GEUS_AWS_position)):
                    if calendar.isleap(GEUS_AWS_position.year[m]):
                        GEUS_AWS_position['n_days']=366
                GEUS_AWS_position['jdy']=GEUS_AWS_position['year']+GEUS_AWS_position['doy']/GEUS_AWS_position['n_days']
       
                plt.plot(GEUS_AWS_position.jdy,GEUS_AWS_position.elev,'s', fillstyle='none',markersize=ms,c='r',
                         label="GEUS GC-Net GPS: %.0f"%np.mean(GEUS_AWS_position.elev)+"±%.0f"%np.std(GEUS_AWS_position.elev)+' m')
                suffix='_has_GEUS_AWS_GPS'
                
                out2=pd.DataFrame({'jdy':GEUS_AWS_position['jdy'],
                                  'elev':GEUS_AWS_position['elev'],
                                  })
                cat_flag=1

            # vals=['t2m']
            # for val in vals:
            #     out_Arctic[val] = out_Arctic[val].map(lambda x: '%.2f' % x)
        
        # if cat_flag:
        #     out_cat = out.append(out2, ignore_index=True)
        # else:
        #     out_cat=out.copy()

        # if len(out_cat)>1:
        #     x=out_cat.jdy.values
        #     y=out_cat.elev.values
        #     v=~np.isnan(y)
        #     y=y[v]
        #     x=x[v]
        #     b, m = polyfit(x, y, 1)
        #     xx=[x[0],x[-1]]
        #     yy=[xx[0]*m+b,xx[1]*m+b]
        #     dy=yy[1]-yy[0]
        #     ny=xx[1]-xx[0]
        #     sign=''
        #     if dy>0:sign='+'
        #     plt.plot(xx,yy,'--',c='grey',
        #              label=f'ATM and GEUS AWS fit\nfrom {"%.1f"%yy[0]}m in {"%.0f"%xx[0]} to {"%.1f"%yy[1]}m in {"%.0f"%int(xx[1])}\n= {sign}{"%.1f"%dy} m elevation change over {"%.0f"%ny} years')

        # out_cat.to_csv(f'./ATM/output/merged_ATM_AWS/{nicknames[k]}.csv',index=None)
        plt.hlines(meta.elev.values[k],year0,year1,linestyle='--',color='grey',label="Table 4 Vandecrux er al 2023: %.0f"%meta.elev.values[k]+' m')

        if nicknames[k]=='JAR':
            GPS_2010=pd.read_csv('./data/GNSS data/2009 JAR1/gnss_jar_all.csv')
            print(GPS_2010.columns)
            geoid_offset=-meta[geoid].values[k]
            GPS_2010['elev']=GPS_2010['ellipsoidal_height_m']+geoid_offset
            GPS_2010['jy']=GPS_2010.year+GPS_2010.day_of_year/365
            mean_GPS_2010=np.mean(GPS_2010['elev'])
            std_GPS_2010=np.std(GPS_2010['elev'])
            plt.plot(GPS_2010['jy'],GPS_2010['elev'],'*',c='r',
                     label="GPS 2010, EGM96: %.0f"%mean_GPS_2010+"±%.0f"%std_GPS_2010+' m',
                     zorder=30)

        if nicknames[k]=='SUM':
            value=3210.442 ; uncert=np.nan
            plt.plot([2011.5,2012.5],[value,value],'-*',c='r',
                     label="Greenland_GNSS_2011_2012.shp %.0f"%value+"±%.0f"%uncert+' m',
                     zorder=30)
            
        if show_approximation_elev_history:
            dy=np.min(elev_approximation)-np.max(elev_approximation)
            dx=np.max(time_elev_approximation)-np.min(time_elev_approximation)
            if elev_approximation[-1]>elev_approximation[0]:dy=-dy
            print('elevation change %.0f'%dy+' m')
            print('N years %.0f'%dx+' m')
            dhdt=dy/dx
            print('linear dhdt %.0f'%dhdt+' m')
            plt.plot(time_elev_approximation,elev_approximation,'-s',c='m',
                     label='approximation of elevation change:\n%.0f'%dy+' m, %.1f'%dhdt+' m y$^{-1}$ over %.1f'%dx+' years',zorder=30)
    
        plt.legend(fontsize=8)
        
        if constrain_elev:
            plt.ylim(elev0,elev1)
        
        # ----------- annotation
        xx0=-0.07 ; yy0=-0.11
        mult=0.8
        co=0.4
        cwd=__file__  #cwd=os.getcwd()
        plt.text(xx0, yy0, cwd.replace('/Users/jason/Dropbox/AWS/GCNET/GCNet_positions','.'),
                fontsize=font_size*mult,color=[co,co,co],rotation=0,
                transform=ax.transAxes,zorder=20,ha='left') # 
        
        ly='p'
        if ly =='x':plt.show()
        if ly =='p':
            plt.savefig(f'./figs/GC-Net_AWS_barometer_elevation_change_reconstructions/{nicknames[k]}_{geoid}.png', bbox_inches='tight', dpi=120)
            # if suffix!='':
            #     plt.savefig(f'./ATM/Figs/has_GEUS_AWS_GPS/{nicknames[k]}_{geoid}_{suffix}.png', bbox_inches='tight', dpi=150)

#%%
wo=0
meta.columns
outx=meta.copy()

if wo:
    outx = outx.rename({'elev': 'elev_assumed_earlier'}, axis=1)
    
    vals=['elev_mean_from_altimetry','elev_linear_intercept','elev_change_linear','elev_change_n_years']
    for val in vals:
        outx[val] = outx[val].map(lambda x: '%.1f' % x)

    vals=['elev_linear_slope']
    for val in vals:
        outx[val] = outx[val].map(lambda x: '%.4f' % x)

    vals=['elev_fit_t0',
                  'elev_fit_t1']
    for val in vals:
        outx[val] = outx[val].map(lambda x: '%.2f' % x)

    header = ['ID', 'name',             'nickname',
'lat', 'lon', 'elev_assumed_earlier', 
              'elev_mean_from_altimetry',
              'elev_linear_slope',
              'elev_linear_intercept',
              'N_altimetry_measurements',
              'elev_fit_t0',
              'elev_fit_t1',
              'elev_change_linear',
              'elev_change_n_years',
              'Date of installation',
           'Last valid time', 'Length of record (years)', 
            # 'k',
            ] 
    outx.to_csv('./ATM/output/GC-Net_elevations_solely_from_ATM_fit.csv', columns = header,index=None)
