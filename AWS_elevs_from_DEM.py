#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Adrien Wehrlé, University of Zurich, Switzerland

"""
import rasterio
import glob
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import matplotlib.animation as animation
import matplotlib as mpl

# prepare paths

if os.getlogin() == "jason":
    base_path = "/Users/jason/Dropbox/rain_optics_SICE_AWS/"

os.chdir(base_path)

#%%
# load geodata

elev = rasterio.open("./metatiffs/elev_1km_1487x2687.tif").read(1)
lat = rasterio.open("./metatiffs/lat_1km_1487x2687.tif").read(1)
lon = rasterio.open("./metatiffs/lon_1km_1487x2687.tif").read(1)
lon[lon==0]=np.nan
lat[lat==0]=np.nan

mask = rasterio.open("./metatiffs/mask_1km_1487x2687.tif").read(1)
# mask[mask>2]=2 # include ice shelves as standard ice mass value of 2
# load profile to further save outputs
profile = rasterio.open("./metatiffs/mask_1km_1487x2687.tif").profile
# plt.hist(mask)
mask[mask == 3] = 2
# mask[mask > 2] = 1
# mask[((mask > 2)&(mask <=200))] = 1
# mask[mask==0]=0

land = np.zeros((2687, 1487))
land = land.astype(np.float64)

ice = np.zeros((2687, 1487))
ice = land.astype(np.float64)

ice[mask == 0] = np.nan
ice[ice==0]=1
tot_ice_area=np.sum(ice==1)

plt.imshow(ice)
plt.colorbar()

#%%
land[mask == 1] = 1
land[mask == 0] = np.nan
land[mask > 1] = 0.5

plt.imshow(land)
plt.colorbar()

# with rasterio.open(
#     "/Users/jason/Dropbox/rain_optics_SICE_AWS/SICE_retrievals/output_tifs/land.tif",
#     "w",
#     **profile
# ) as dst:
#     dst.write(land, 1)

# %% map the good stuff

meta=pd.read_csv('/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/output/GC-Net_info_incl_1999.csv')
print(meta.columns)
print(meta)

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


n=len(meta)

meta['elev_from_dem']=np.nan

for st in range(n):
# for st,site in enumerate(sites):
    # if st>=0:
    # if site=='tunu_n':
    # if site=='summit':
    print(st,meta.name[st])
    if meta.name[st]!='null':
    # if meta.name[st]=='DY2':
        dist=haversine_np(meta.lon[st],meta.lat[st], lon, lat)
        # plt.imshow(dist)
        # plt.colorbar()
        dist[np.isnan(dist)]=1000
        v=np.where(dist==np.min(dist))
        print(lat[v],lon[v],elev[v])
        meta['elev_from_dem'][st]=int(elev[v][0])

meta.to_csv('/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/output/GCNet_sites_with_elev_from_DEM.csv')
meta.to_csv('/Users/jason/Dropbox/AWS/CARRA-TU_GEUS/metadata/GCNet_sites_with_elev_from_DEM.csv')
#%%
passes=['ASC','DSC']
times=['13','01']
for j,passx in enumerate(passes):
    # files = sorted(glob.glob("/Users/jason/0_dat/AMSR2/*"+passx+"*.tif"))
    files = sorted(glob.glob("/Users/jason/0_dat/AMSR2/melt_geotif/*"+passx+"*.tif"))
    print(files)
    print(len(files))
    x=np.array(files)
    print(passx,files)
    
    # 13 UTC ascending and 01 UTC descending passes
    font_size=16
    for i,file in enumerate(files):
        # if i>=0:
        if i==0:
            img = rasterio.open(file).read(1)
            img = img.astype(np.float64)
            
            
            img[mask == 0] = np.nan
            img[land==1] = np.nan
    
            date=x[i].split('/')[-1][9:19]+' '+times[j]+'UTC'
            print(date)
            
            plt.close()
            plt.clf()
            fig, ax = plt.subplots(figsize=(5, 5))
    
            ax.imshow(img,cmap='Reds')
            ax.axis('off')
            
            cc=0
            xx0=0.01 ; yy0=0.0 ; dy2=-0.04
            mult=0.7
            color_code='grey'
            plt.text(xx0, yy0+cc*dy2, 'Picard & Fily 2006',
                    fontsize=font_size*mult,color=color_code,rotation=0,transform=ax.transAxes) ; cc+=1.5
    
            xx0=0.01 ; yy0=1.05 ; dy2=-0.04
            mult=0.7
            color_code='grey'
            plt.text(xx0, yy0+cc*dy2, date,
                    fontsize=font_size*mult,color=color_code,rotation=0,transform=ax.transAxes) ; cc+=1.5
    
            # plt.title(date)
    
            ly = "p"
        
            if ly == "p":
                plt.savefig('./AMSR2/figs/'+date+'_v2.png', bbox_inches="tight")

    #%%

    # print(df)
    # # %%

    # dalbedo_prof_perc = ((np.array(albedo_prof_poevent) - np.array(albedo_prof_prevent))\
    #     / np.array(albedo_prof_prevent)) * 100

    # plt.figure(figsize=(13, 9))
    # plt.plot(elev_bins[:-1], dalbedo_prof_perc, color='gray')
    # plt.scatter(elev_bins[:-1], dalbedo_prof_perc, color='gray', alpha=0.5)
    # plt.xlim([300, 2150])
    # plt.grid(alpha=0.5)
    # plt.tick_params(labelsize=18)
    # plt.xlabel('Elevation (meters)', fontsize=20)
    # plt.ylabel('Albedo decrease (%)', fontsize=20)

    # plt.axvline(elev_bins[:-1][np.array(albedo_prof_prevent) >= 0.565][0], color=color0,
    #             LineStyle='--', label='Pre-event snowline elevation')

    # plt.axvline(elev_bins[:-1][np.array(albedo_prof_poevent) >= 0.565][0], color=color1,
    #             LineStyle='--', label='Post-event snowline elevation')
    # plt.legend(fontsize=15)
    # plt.title('Albedo decrease/Elevation profiles in CW Greenland \n before and after extreme 2021-08 rainfall event',
    #           fontsize=23)

    # # plt.savefig('.figures/082021_albedo_elev_profile_CW_perc_decrease_perc.png',
    # #             bbox_inches='tight')

    # # %%

    # dalbedo_prof = np.array(albedo_prof_poevent) - np.array(albedo_prof_prevent)

    # plt.figure(figsize=(13, 9))
    # plt.plot(elev_bins[:-1], dalbedo_prof, color='gray')
    # plt.scatter(elev_bins[:-1], dalbedo_prof, color='gray', alpha=0.5)
    # plt.xlim([300, 2150])
    # plt.grid(alpha=0.5)
    # plt.tick_params(labelsize=18)
    # plt.xlabel('Elevation (meters)', fontsize=20)
    # plt.ylabel('Albedo decrease (unitless)', fontsize=20)

    # plt.axvline(elev_bins[:-1][np.array(albedo_prof_prevent) >= 0.565][0], color=color0,
    #             LineStyle='--', label='Pre-event snowline elevation')

    # plt.axvline(elev_bins[:-1][np.array(albedo_prof_poevent) >= 0.565][0], color=color1,
    #             LineStyle='--', label='Post-event snowline elevation')
    # plt.legend(fontsize=15)
    # plt.title('Albedo decrease/Elevation profiles in CW Greenland \n before and after extreme 2021-08 rainfall event',
    #           fontsize=23)

    # # plt.savefig('./figures/082021_albedo_elev_profile_CW_perc_decrease.png',
    # #             bbox_inches='tight')
