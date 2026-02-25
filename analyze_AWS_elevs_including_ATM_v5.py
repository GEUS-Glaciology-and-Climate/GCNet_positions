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

the code outputs position_interpolated_with_elev.csv files used by ./GC-Net_1995-2020_w_hourly_positions/src/merge_positions_incl_elev_with_met_data.py

"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pathlib import Path
import calendar


BASE_PATH = Path("/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/") if os.getlogin() == "jason" else Path(".")

meta = (
    pd.read_csv(BASE_PATH / "meta/GC-Net_location_essd.csv")
    .rename(
        columns={
            "Name": "name",
            "Latitude (°N)": "lat",
            "Longitude (°E)": "lon",
            "Elevation (wgs84 m)": "elev",
        }
    )
)
nicknames= [ "SWC", "CP1", "CP2", "JAR", "JR2", "JR3",
     "NAU", "GIT", "HUM", "SUM", "TUN", "DY2",
     "SDL", "SDM", "NAE", "NSE", "NGRP", "NEM",
     "EGP", "KAR", "KUL", "AUR", "PTG",
     "PET",  ]
nickname_map = dict(
    zip(
        [
            "Swiss Camp", "Crawford Point 1", "CP2", "JAR1", "JAR2", "JAR3",
            "NASA-U", "GITS", "Humboldt", "Summit", "Tunu-N", "DYE-2",
            "Saddle", "South Dome", "NASA-E", "NASA-SE", "NGRIP", "NEEM",
            "EastGRIP", "KAR", "KULU", "Aurora", "Petermann Glacier",
            "Petermann ELA",
        ],
       nicknames,
    )
)

meta["nickname"] = meta["name"].map(nickname_map).fillna("")

drop_names = {
    "SMS-PET","SMS1","SMS2","SMS3","SMS4","SMS5","LAR1","LAR2","LAR3",
    "Swiss Camp 10m","Roof_GEUS","LYN_T","LYN_L","CEN1","CEN2",
    "KPC_Lv3","KPC_Uv3","THU_L2","WEG_B","ZAK_Uv3",
}
meta = meta.loc[~meta["name"].isin(drop_names)].reset_index(drop=True)

new_cols = [
    "elev_linear_slope",
    "elev_linear_intercept",
    "elev_fit_t0",
    "elev_fit_t1",
    "elev_change_linear",
    "elev_change_n_years",
    "elev_mean_from_altimetry",
    "k",
    "N_altimetry_measurements",
]
meta[new_cols] = np.nan

print(meta["name"])
print(meta.columns)

import numpy as np
import pandas as pd
import rasterio

GEOIDS = ["egm96_15"]
EGM_TIF = "meta/us_nga_egm96_15.tif"

for geoid in GEOIDS:
    with rasterio.open(EGM_TIF) as ds:
        pts = list(zip(meta["lon"].astype(float).to_numpy(), meta["lat"].astype(float).to_numpy()))
        meta[geoid] = np.fromiter((v[0] for v in ds.sample(pts)), dtype="float64", count=len(pts))

df = (
    pd.read_excel("./meta/GC-Net historical positions.xlsx")
    .query("date not in ['doc_2000','ref','WMO-DMI_2012']")
    .loc[lambda d: d["date"] != 35596]
    .loc[lambda d: d["elev"].ne("-")]
    .assign(
        name_long=lambda d: d["name_long"].replace({"Crawford Pt. 1": "Crawford Point 1"}),
        date=lambda d: pd.to_datetime(d["date"], errors="coerce"),
    )
    .dropna(subset=["date"])
)

df = df.assign(
    year=df["date"].dt.year,
    doy=df["date"].dt.dayofyear,
    jy=df["date"].dt.year + df["date"].dt.dayofyear / 366.0,
)

vandecrux_position_compilation = df.copy()

fn = "./data/GNSS data/2023/Data_overview_2023.xlsx"
GNSS_2023 = (
    pd.read_excel(fn, sheet_name="Calculated_positions", skiprows=1)
    .rename(columns={"Station name": "name", "Ortho H": "elev"})
    .assign(date=pd.to_datetime("2023-06-15"))
    .assign(
        year=lambda d: d["date"].dt.year,
        doy=lambda d: d["date"].dt.dayofyear,
        jy=lambda d: d["date"].dt.year + d["date"].dt.dayofyear / 366.0,
    )
)

n_AWS = len(meta)
iyear, fyear = 1995, 1998
years = np.arange(iyear, fyear + 1).astype(str)
n_years = len(years)

# --- plot style ---
th = 1
ms = 12
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.facecolor": "w",
        "axes.edgecolor": "k",
        "axes.grid": True,
        "grid.alpha": 1,
        "grid.color": "#cccccc",
        "grid.linewidth": th,
        "axes.linewidth": th,
        "legend.facecolor": "w",
        "legend.framealpha": 0.8,
        "mathtext.default": "regular",
        "figure.figsize": (5, 4),
    }
)

sites = [
    "Swiss Camp", "Crawford Point 1", "CP2", "JAR1", "JAR2", "JAR3",
    "NASA-U", "GITS", "Humboldt", "Summit", "Tunu-N", "DYE-2",
    "Saddle", "South Dome", "NASA-E", "NASA-SE", "NGRIP", "NEEM",
    "EastGRIP", "KAR", "KULU", "Aurora", "Petermann Glacier", "Petermann ELA",
]
sites2 = [
    "Swiss Camp", "Crawford Pt. 1", "CP2", "JAR1", "JAR2", "JAR3",
    "NASA-U", "GITS", "Humboldt", "Summit", "Tunu-N", "DYE-2",
    "Saddle", "South Dome", "NASA-E", "NASA-SE", "NGRIP", "NEEM",
    "East Grip", "KAR", "KULU", "Aurora", "Petermann Glacier", "Petermann ELA",
]

SPECIAL_FLAGS = {"SWC": {"show_stober": 1, "show_khan": 1}, "JAR": {"show_stober": 1, "show_khan": 1}}
NO_INTERP = {"GIT", "SUM", "SDL", "SDM", "NGRP", "KAR", "JR2","JR3","KUL","AUR","PTG","PET"}

MIN_TOL = {
    "HUM": 8, "NSE": 4, "NGRP": 4, "NEM": 4, "EGP": 4, "KAR": 4,
    "KUL": 2, "AUR": 5, "PTG": 5, "PET": 3,
}

ELEV_RULES = {
    "SWC": ([1990.0, 1990.5, 2005, 2015, 2020, 2023.5], [1155.5, 1155.5, 1135, 1122, 1119.2, 1119.2]),
    "CP1": ([1995, 2023.5], ["FIX_1960", "GNSS"]),
    "JAR": ([1996.45, 2010, 2015, 2023.5], [932.3, 919.5, 909, 906]),
    "JR2": ([1999, 2020], [523, 484]),
    "JR3": ([2000, 2020], [300, 220]),
    "NAU": ([1995, 2023.5], [2338, "GNSS"]),
    "GIT": ([1995, 2022], [1872, 1871]),
    "HUM": ([1995, 2023.5], [1974, "GNSS"]),
    "SUM": ([1996, 2023], [3205, 3205]),
    "TUN": ([1996, 2023], [2074, 2078]),
    "DY2": ([1996, 2023.5], [2115, "GNSS"]),
    "SDL": ([1997, 2023.5], [2456, "GNSS"]),
    "SDM": ([1997, 2023.5], [2878, "GNSS"]),
    "NAE": ([1997, 2023.5], [2622.7, "GNSS"]),
    "NSE": ([1998, 2023.5], [2372, "GNSS"]),
    "NGRP": ([1997, 2021], [2919.3, 2919.3]),
    "NEM": ([2006, 2023.5], [2451, 2450.2]),
    "EGP": ([2014, 2023.5], [2663.5, 2663.2]),
    "KAR": ([1999, 2002], [2552.5, 2552.5]),
    "KUL": ([1999, 2001], [900, 900]),
    "AUR": ([1999, 2001], [1748, 1748]),
    "PTG": ([2002, 2007], [33, 28]),
    "PET": ([2002, 2021], [949, 940]),
}

def elev_v_time_interpolation(fn, site_out, time_elev_approximation, elev_approximation):
    pos = pd.read_csv(fn, parse_dates=["date"])
    year = pos["date"].dt.year
    doy = pos["date"].dt.dayofyear
    ndays = np.where(year.map(calendar.isleap).to_numpy(), 366.0, 365.0)
    jdy = year.to_numpy() + doy.to_numpy() / ndays

    t = np.asarray(time_elev_approximation, dtype="float64")
    z = np.asarray(elev_approximation, dtype="float64")
    pos["elev"] = np.interp(jdy, t, z, left=np.nan, right=np.nan)

    pos["elev"] = pos["elev"].round(1)
    pos.to_csv(f"./output/{site_out}_position_interpolated_with_elev.csv", columns=["date", "lon", "lat", "elev"], index=False)


def load_ATM(nick, geoid_offset):
    ATM = pd.read_csv(f"./ATM/output/{nick}.csv")
    ATM = add_jdy(ATM, "date")
    elev = ATM["elev_ATM_m"] + geoid_offset
    dist = ATM["distance_km"].to_numpy()
    x1 = ATM["slope_S2N_deg"].to_numpy() * dist * 1000
    x2 = ATM["slope_W2E_deg"].to_numpy() * dist * 1000
    yx = ATM["elev_ATM_m"].to_numpy() + geoid_offset - x1 - x2
    return ATM, elev.to_numpy(), dist, yx

BAD_ATM_FILTERS = {
    "CP1": lambda t, y: (t > 2005) & (y > 1960),
    "CP2": lambda t, y: (t > 2005) & (y > 1960),
    "JR2": lambda t, y: (t > 2010) & (y > 530),
    "JR3": lambda t, y: (t > 2005) & (y < 220),
}

def add_jdy(df, date_col="date"):
    d = pd.to_datetime(df[date_col])
    year = d.dt.year
    doy = d.dt.dayofyear
    ndays = np.where(year.map(calendar.isleap).to_numpy(), 366.0, 365.0)
    df = df.copy()
    df[date_col] = d
    df["jdy"] = year.to_numpy() + doy.to_numpy() / ndays
    return df

def read_stober_xlsx(fn, drop_head=0, drop_tail=0):
    st = pd.read_excel(
        fn,
        skiprows=7,
        names=["latdeg","latmin","latsec","NS","lat","emt","londeg","lonmin","lonsec","EW","lon","elevation"],
        index_col=0,
    )
    if drop_head:
        st = st.iloc[drop_head:]
    if drop_tail:
        st = st.iloc[:-drop_tail]
    st = st.reset_index(names="raw").astype({"raw": "string"})
    return st

def parse_ddmmyy(s):
    return pd.to_datetime(s, format="%d%m%y", errors="coerce")

def load_stober_JAR(fn, geoid_offset, ST2_offset=35):
    st = read_stober_xlsx(fn, drop_head=65)
    stake = st["raw"].str.slice(0, 5)
    date = parse_ddmmyy(st["raw"].str.slice(6))
    m = (stake == "ST201") & date.notna()
    out = pd.DataFrame(
        {
            "date": date[m],
            "elev": st.loc[m, "elevation"].to_numpy() - 1.08 + geoid_offset - ST2_offset,
        }
    )
    return add_jdy(out)

def load_stober_SWC(fn, geoid_offset, peg_ids=("106","120","121"), geoid_override=-25):
    st = read_stober_xlsx(fn, drop_tail=32)
    pegel = st["raw"].str.slice(0, 3)
    date = parse_ddmmyy(st["raw"].str.slice(4))
    m = date.notna() & (st["raw"].str.slice(4) != "ex010915")
    st = st.loc[m].assign(pegel=pegel[m].to_numpy(), date=date[m].to_numpy())

    go = geoid_override if geoid_override is not None else geoid_offset
    elev = st["elevation"].to_numpy() - 1.08 + go

    out = (
        pd.DataFrame({"date": st["date"], "pegel": st["pegel"], "elev": elev})
        .query("pegel in @peg_ids")
        .pivot_table(index="date", columns="pegel", values="elev", aggfunc="mean")
        .rename(columns={"106": "elev0", "120": "elev1", "121": "elev2"})
        .reset_index()
    )
    return add_jdy(out)

stober_JAR = load_stober_JAR(
        "./meta/Stober_et_al_2023/SWC+ST2-Alle Pegel_Koordinaten Geodätisch.xlsx",
        geoid_offset=-meta.loc[meta.nickname=='JAR',geoid].iloc[0],
        ST2_offset=35,
    )
stober_SWC = load_stober_SWC(
        "./meta/Stober_et_al_2023/SWC+ST2-Alle Pegel_Koordinaten Geodätisch.xlsx",
        geoid_offset=-meta.loc[meta.nickname=='SWC',geoid].iloc[0],
        geoid_override=-25,
    )

df_meta_new = pd.read_csv('meta/GC-Net_location_essd.csv')
# %% --- per-site ---
df_list = []
for nick in ELEV_RULES.keys():
    print(nick)
    k = np.argwhere(np.array(nicknames)==nick)[0][0]
    site_out = sites2[k]

    show_khan = SPECIAL_FLAGS.get(nick, {}).get("show_khan", 0)
    show_stober = SPECIAL_FLAGS.get(nick, {}).get("show_stober", 0)

    min_tolerated_dist = MIN_TOL.get(nick, 1)

    GNSS_2023_elev_pick = np.nan
    GNSS_2023_decimal_year = np.nan
    v = (GNSS_2023["name"] == nick)
    if v.any():
        GNSS_2023_elev_pick = float(GNSS_2023.loc[v, "elev"].iloc[0])
        GNSS_2023_decimal_year = float(GNSS_2023.loc[v, "jy"].iloc[0])

    time_elev_approximation, elev_approximation = ELEV_RULES[nick]
    elev_approximation = [
        (1960.0 if v == "FIX_1960" else GNSS_2023_elev_pick if v == "GNSS" else float(v))
        for v in elev_approximation
    ]

    geoid_offset = -float(meta[geoid].values[k])

    fig, ax = plt.subplots(figsize=(9, 7))

    ATM, elevs, dist, yx = load_ATM(nick, geoid_offset)
    v = np.where(dist < min_tolerated_dist)[0]

    ATM_dates = ATM["jdy"].to_numpy()[v]
    elevs_tol = elevs[v].astype(float)

    bad = BAD_ATM_FILTERS.get(nick)
    if bad is not None:
        elevs_tol[bad(ATM_dates, elevs_tol)] = np.nan

    lab = f"ATM within {min_tolerated_dist:.1f} km: {np.nanmean(elevs_tol):.1f}±{np.nanstd(elevs_tol):.1f} m, using {geoid.split('_')[0].upper()}"
    ax.plot(ATM_dates, elevs_tol, "o", fillstyle="none", markersize=ms, mew=2, label=lab)

    # --- Khan/Stober overlays ---
    if nick == "SWC":
        lab = f"GPS survey c/o M. Stober, Jan 2024, stakes 106: {np.nanmean(stober_SWC['elev0']):.1f}±{np.nanstd(stober_SWC['elev0']):.1f} m"
        ax.plot(stober_SWC["jdy"], stober_SWC["elev0"], "s", color="k", label=lab)
    elif nick == "JAR":
        lab = f"Stober ST201 incl. 35 m offset: {np.nanmean(stober_JAR['elev']):.1f}±{np.nanstd(stober_JAR['elev']):.1f} m"
        ax.plot(stober_JAR["jdy"], stober_JAR["elev"], "s", color="k", zorder=20, label=lab)

    if nick in ["SWC", "JAR"]:
        khan_fn = "./Khan/SwissCamp.txt" if nick == "SWC" else "./Khan/JAR.txt"
        khan = pd.read_csv(khan_fn, skiprows=2, sep=r"\s+", names=["jy", "elev"])
        os_khan = 1135.7 if nick == "SWC" else 927.0
        ax.plot(khan["jy"], khan["elev"] + os_khan, "-", linewidth=th * 4, color="orange", zorder=20,
                label="satellite altimetry c/o S.A. Khan, Jan. 2024\noffset to GEUS GPS data 2022 to 2023")

    # --- GNSS 2023 ---
    if not np.isnan(GNSS_2023_elev_pick):
        ax.plot(GNSS_2023_decimal_year, GNSS_2023_elev_pick, "o", fillstyle="none", color="c",
                markersize=ms, mew=2, zorder=20, label=f"Jakobsen et al. GEUS 2023 GNSS: {GNSS_2023_elev_pick:.1f}±0.1 m")

    # --- historical positions ---
    vh = vandecrux_position_compilation["name_long"].eq(sites[k])
    y3 = vandecrux_position_compilation.loc[vh, "elev"].to_numpy(dtype="float64")
    if np.isfinite(np.nanstd(y3)):
        ax.plot(vandecrux_position_compilation.loc[vh, "jy"], y3, "s", fillstyle="none", markersize=ms / 2, color="b",
                label=f"GC-Net historical positions.xlsx: {np.nanmean(y3):.0f}±{np.nanstd(y3):.0f} m")

    # --- GEUS AWS GPS (monthly) ---
    fn = Path(f"./output/Jason/{nick}_positions_monthly.csv")
    if fn.is_file():
        gps = pd.read_csv(fn)
        gps.loc[gps["elev"] == 0, "elev"] = np.nan
        gps["elev"] = gps["elev"] - 1.5
        gps["date"] = pd.to_datetime(dict(year=gps["year"], month=gps["month"], day=15))
        gps = add_jdy(gps, "date")
        ax.plot(gps["jdy"], gps["elev"], "s", fillstyle="none", markersize=ms, color="r",
                label=f"GEUS GC-Net GPS: {np.nanmean(gps['elev']):.0f}±{np.nanstd(gps['elev']):.0f} m")

    # --- meta line ---
    year0, year1 = 1989, 2024
    ax.hlines(meta.elev.values[k], year0, year1, linestyle="--", color="grey",
              label=f"Table 4 Vandecrux et al 2023: {meta.elev.values[k]:.0f} m")

    # --- JAR 2010 GNSS ---
    if nick == "JAR":
        GPS_2010 = pd.read_csv("./data/GNSS data/2009 JAR1/gnss_jar_all.csv")
        GPS_2010["elev"] = GPS_2010["ellipsoidal_height_m"].to_numpy() + geoid_offset
        GPS_2010["jy"] = GPS_2010["year"].to_numpy() + GPS_2010["day_of_year"].to_numpy() / 365.0
        ax.plot(GPS_2010["jy"], GPS_2010["elev"], "*", color="r", zorder=30,
                label=f"GPS 2010, EGM96: {GPS_2010['elev'].mean():.0f}±{GPS_2010['elev'].std():.0f} m")

    if nicknames[k]=='SUM':
        value=3210.442 ; uncert=np.nan
        plt.plot([2011.5,2012.5],[value,value],'-*',c='r',
                 label="Greenland_GNSS_2011_2012.shp %.0f"%value+"±%.0f"%uncert+' m',
                 zorder=30)

    dy=np.min(elev_approximation)-np.max(elev_approximation)
    dx=np.max(time_elev_approximation)-np.min(time_elev_approximation)
    if elev_approximation[-1]>elev_approximation[0]:dy=-dy
    print('elevation change %.0f'%dy+' m')
    print('N years %.0f'%dx+' m')
    dhdt=dy/dx
    print('linear dhdt %.0f'%dhdt+' m')
    plt.plot(time_elev_approximation,elev_approximation,'-s',c='m',
             label='approximation of elevation change:\n%.0f'%dy+' m, %.1f'%dhdt+' m y$^{-1}$ over %.1f'%dx+' years',zorder=30)
    df_meta_new.loc[df_meta_new.Name == sites[k], 'Elevation (wgs84 m)'] = np.round(np.mean(elev_approximation))

    ax.hlines(df_meta_new.loc[df_meta_new.Name == sites[k], 'Elevation (wgs84 m)'] ,
              year0, year1, linestyle="--", color="tab:green",
              label=f"New elevation in GC-Net Level 1 metadata: {meta.elev.values[k]:.0f} m")
    df_site = pd.DataFrame()
    df_site['year'] = time_elev_approximation
    df_site['altitude'] = elev_approximation
    df_site['site'] = sites[k]

    df_list.append(df_site)

    # --- titles / axes / save ---
    ax.set_title(f"{sites[k]} a.k.a. {nick}")
    ax.set_ylabel("AWS barometer elevation above mean sea level, m")
    ax.set_xlim(year0, year1)
    ax.legend(fontsize=8)

    plt.savefig(f"./figs/GC-Net_AWS_barometer_elevation_change_reconstructions/{nick}_{geoid}.png",
                bbox_inches="tight", dpi=120)

pd.concat(df_list).to_csv('output/GC-Net_elevation_tie_points.csv', index=None)
df_meta_new.to_csv('output/GC-Net_locations.csv', index=None)
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
