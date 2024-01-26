#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 10:22:26 2023

@author: jason, jeb@geus.dk

gather data using

https://nsidc.org/icebridge/portal/map
via
https://nsidc.org/data/ilatm1b/versions/2

select checkboxes for granules for either the old ATM (BLATM2) or the Icebridge ATM data (ILATM2),
I have archived the .zip files, takes quite some time to gather

"""

from glob import glob
import pandas as pd
import numpy as np
import os

# ## change to your system's login name to change dir for local work
if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/AWS/GCNET/GCNet_positions.stash/'

os.chdir(base_path)

def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.    

    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km

def lon360_to_lon180(lon360):

    # reduce the angle  
    lon180 =  lon360 % 360 
    
    # force it to be the positive remainder, so that 0 <= angle < 360  
    lon180 = (lon180 + 360) % 360;  
    
    # force into the minimum absolute value residue class, so that -180 < angle <= 180  
    lon180[lon180 > 180] -= 360
    
    return lon180


# read Greenland AWS locations
meta = pd.read_csv('/Users/jason/Dropbox/AWS/GCNET/GC-Net-level-1-data-processing/L1/GC-Net_location.csv')
meta = meta.rename({'Name': 'name'}, axis=1)
meta = meta.rename({'Latitude (°N)': 'lat'}, axis=1)
meta = meta.rename({'Longitude (°E)': 'lon'}, axis=1)
meta = meta.rename({'Elevation (wgs84 m)': 'elev'}, axis=1)
meta = meta.rename({'Crawford Point 1': 'Crawford Pt. 1'}, axis=1)
meta=meta.replace({'Crawford Point 1': 'Crawford Pt. 1'})


names=['SMS-PET','SMS1','SMS2','SMS3','SMS4','SMS5','LAR1','LAR2','LAR3','Swiss Camp 10m','Roof_GEUS','LYN_T','LYN_L','CEN1','CEN2','KPC_Lv3','KPC_Uv3','THU_L2','WEG_B','ZAK_Uv3']
for name in names:
    meta.drop(meta[meta.name==name].index, inplace=True) # drop original time column

print(meta.name)
# print(meta.columns)


#%%
n_AWS=len(meta)

# meta=meta.sort_values(by='name', ascending=True)
# print(meta.columns)

min_tolerated_dist=20

iyear=1995 ; fyear=1995
iyear=1996 ; fyear=1999 ; year_prefix=1900
iyear=2001 ; fyear=iyear ; 
n_years=fyear-iyear+1
years=np.arange(iyear,fyear+1).astype(str)

years=['1995',
'1996',
'1997',
'1998',
'1999',
'2001',
'2002',
'2003',
'2005',
'2006',
'2008',
]


years=[
       
# '2003',
'2005',
'2006',
'2008',
]

# years=['2010',
#        '2011',
#        '2012',
#        '2013',
#        '2014',
#        '2015',
#        '2016',
#        '2017',
#        '2019',
#        ]

# years=[
#         '2003',
#         ]

n_years=len(years)

for year in years:
    print(year)
    if int(year)>=2000:
        year_prefix=2000
    else:
        year_prefix=1900
    
    path='/Users/jason/0_dat/ATM/'+year+'/'
    
    files = glob(path+"*")
    
    print(files)
    
    n_files=len(files)
    
    for i,file in enumerate(files):

        print(year,i,n_files-i)
        if int(year)<2009:
            fn = glob(file+"/*seg*")[0]
            temp=fn.split('/')[-1]
            temp2=temp.split('_')[1]
            yearx=int(temp2[0:2])
        else:
            fn = glob(file+"/ILATM2*.csv")[0]
            temp=fn.split('/')[-1]
            temp2=temp.split('_')[1]
            yearx=int(temp2[0:4])

        if int(year)<2009:
            if yearx>10:
                yearx+=1900
            else:
                yearx+=year_prefix

            month=temp2[2:4]
            day=temp2[4:6]
        else:
            month=temp2[4:6]
            day=temp2[6:8]

        if yearx!=int(year):
            print(fn)
            print(yearx,int(year))
            os.system('rm -R '+file)

        if yearx==int(year):
            if i>=0:
            # if i<=5:
            # if ((i>=5)&(i<=10)):
            # if i==3:
                # flag=0
                
                elevs=np.zeros(n_AWS)
                dists=np.zeros(n_AWS)+900
                lats=np.zeros(n_AWS)
                lons=np.zeros(n_AWS)
                slope_S2N=np.zeros(n_AWS)
                slope_W2E=np.zeros(n_AWS)
                yearsx=np.zeros(n_AWS)
                monthsx=np.zeros(n_AWS)
                daysx=np.zeros(n_AWS)
                if int(year)<2009:
                    ATM_raw=pd.read_csv(fn,delim_whitespace=True,names=['time','cenlat','cenlon','z','S2N','W2E','rms','N','Nbad','dist','block'])
                else:
                    # print('hi')
                    ATM_raw=pd.read_csv(fn,names=['time','cenlat','cenlon','z','S2N','W2E','rms','N','Nbad','dist','block'],skiprows=11)
                # print(fn)
                ATM_raw['cenlon180']=lon360_to_lon180(ATM_raw.cenlon)

                # print(year,n_files,i,n_files-i,fn,yearx,month,day)
            
                
                count_ok=0
                for k,site in enumerate(meta.name):
                # for k in range(1):
                    if k>=0:
                        # print(k,site)
                        if site!='CP2' and site!='JAR2'and site!='JAR3'and site!='GITS'and site!='Summit'and site!='Saddle'and site!='South Dome'and site!='NGRIP'and site!='NEEM'and site!='EastGRIP'and site!='KAR'and site!='KULU'and site!='Aurora'and site!='Petermann Glacier':
                            pos = pd.read_csv(f'/Users/jason/Dropbox/AWS/GCNET/GCNet_positions.stash/output/{site}_position_interpolated.csv')
                            pos.columns
                            
                            pos.date=pd.to_datetime(pos.date).dt.strftime('%Y%m%d')
                            # print(pos.date)
                            
                            v=np.where(pos.date==str(yearx)+day+month)
                            
                            if v[0].size>0:
                                target_lon,target_lat=pos.lon.values[v[0][0]],pos.lat.values[v[0][0]]
                            else:
                                target_lon,target_lat=meta.lon.values[k],meta.lat.values[k]
                                
                            # print(target_lon,target_lat,meta.lon.values[k],meta.lat.values[k])
                        else:
                            target_lon,target_lat=meta.lon.values[k],meta.lat.values[k]
                        # print(k,site,target_lon,target_lat)

                    

                        min_dist=900
                        elev_ATM=0
                        # dist=haversine_np(meta.lon.values[k],meta.lat.values[k], ATM_raw.cenlon180.values, ATM_raw.cenlat.values)
                        dist=haversine_np(target_lon,target_lat, ATM_raw.cenlon180.values, ATM_raw.cenlat.values)
                        min_dist=np.min(dist)
                        v=np.where(dist==np.min(dist))
                        v=v[0]
                        if len(v)>1:v=v[0]
                        # print(np.shape(dist),min_dist)
                        
                        elev_ATM=ATM_raw.z[v]
                        if min_dist<min_tolerated_dist:
                            # print(year,meta.name.values[k])
                            count_ok+=1
                        
                        # if min_dist<dists[k]:
                            elevs[k]=ATM_raw.z[v]
                            dists[k]=min_dist
                            slope_S2N[k]=ATM_raw.S2N[v]
                            slope_W2E[k]=ATM_raw.W2E[v]
                            lats[k]=ATM_raw.cenlat[v]
                            lons[k]=ATM_raw.cenlon180[v]
                            yearsx[k]=yearx
                            monthsx[k]=int(month)
                            daysx[k]=int(day)
                            # sentence=meta.name.values[k],elev_ATM,min_dist
                            # sentence_list.append(sentence)
                            # flag=1
    
                    # if count_ok==0:
                        # print('no hits')
                        # os.system('rm -R '+file)
                    
            out=pd.DataFrame({'site':meta.name,
                          'year':yearsx,
                          'month':monthsx,
                          'day':daysx,
                          'elev_ATM':elevs,
                          'dist':dists,
                          'slope_S2N':slope_S2N,
                          'slope_W2E':slope_W2E,
                          'lat':lats,
                          'lon':lons
                              })
            
            vals=['elev_ATM']
            for val in vals:
                out[val] = out[val].map(lambda x: '%.2f' % x)
                
            vals=['dist']
            for val in vals:
                out[val] = out[val].map(lambda x: '%.3f' % x)

            vals=['lat','lon']
            for val in vals:
                out[val] = out[val].map(lambda x: '%.4f' % x)

            vals=['slope_W2E','slope_S2N']
            for val in vals:
                out[val] = out[val].map(lambda x: '%.4f' % x)
            
            out.to_csv(f'./ATM/output/{year}_{file.split(os.sep)[-1]}.csv')
            # print(out)
                        # plt.scatter(df.cenlat,df.z)
                        # df.time
        