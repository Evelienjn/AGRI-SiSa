#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 09:32:29 2023

@author: evelien

The goal of this script is to calculate spei and spei maps

eraland data
"""
import rpy2

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from tqdm import tqdm
from utils import *
import xarray as xr
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
#%%
r_time_series = robjects.r('ts')
SPEI_package = importr('SPEI')
r_spei_function = robjects.r['spei']
r_spi_function = robjects.r['spi'] 

tstart = '1990-01'
tfinal='2022-12'


data_path_era5land = '/home/evelien/Agri-SiSa/data/ERA5/ERA5land_tp_ev_soilmoisture_mm_1959_2022_NL.nc' 
#select expver=1 (niet erat (prilimilarly), maar era5=validated!)
ds= xr.open_dataset(data_path_era5land).sel(time=(slice(tstart, tfinal)), expver=1)

#%%
#grid settings
longitude_dim =ds.coords['longitude'].values
latitude_dim= ds.coords['latitude'].values
time_dim = ds.coords['time'].values
xx_lon, yy_lat = np.meshgrid(longitude_dim,
                   latitude_dim,
                   )

xdim=ds.coords['longitude'].shape[0]
ydim=ds.coords['latitude'].shape[0]
tdim=ds.coords['time'].shape[0]

spi1_tyx_empty = np.full(shape=[tdim, ydim, xdim],
                 fill_value=np.nan)
spi3_tyx_empty = np.full(shape=[tdim, ydim, xdim],
                 fill_value=np.nan)
spi6_tyx_empty = np.full(shape=[tdim, ydim, xdim],
                 fill_value=np.nan)
spei1_tyx_empty = np.full(shape=[tdim, ydim, xdim],
                 fill_value=np.nan)
spei3_tyx_empty = np.full(shape=[tdim, ydim, xdim],
                 fill_value=np.nan)
spei6_tyx_empty = np.full(shape=[tdim, ydim, xdim],
                 fill_value=np.nan)

#%%
#oefenen met indexen van 1 tijdserie op een coordinaten punt



def pr_at_index(ds, x_ind, y_ind):
    
    ds2=ds.isel(longitude=x_ind,
        latitude= y_ind)
    da = ds2.tp * 1000
    da.attrs={'units': 'mm', 'long_name': 'Total precipitation', 'data_source':'Era5Land_expver1' }
    
    return da


def pev_at_index(ds, x_ind, y_ind):
    ds2= ds.isel(longitude=x_ind,
                 latitude=y_ind
                )
    da= ds2.pev * -1000
    da.attrs={'units': 'mm', 'long_name': 'Potential evaporation', 'data_source':'Era5Land_expver1'}
    return da

#%%

for x in tqdm(np.arange(xdim)):
    for y in np.arange(ydim):
        
        print(longitude_dim[x],' ', latitude_dim[y])
        
        da_pr = pr_at_index(ds=ds, x_ind=x, y_ind=y)
        da_pet = pev_at_index(ds=ds, x_ind=x, y_ind=y)
        #spi_yx = SPI_1d(da_pr=da_pr)
        spei_yx = SPEI_1d(da_pr=da_pr, da_pet=da_pet)
        
        #spi1_tyx_empty[:, y, x] = spi_yx['SPI1'].values
        #spi3_tyx_empty[:, y, x] = spi_yx['SPI3'].values
        #spi6_tyx_empty[:, y, x] = spi_yx['SPI6'].values
        
        spei1_tyx_empty[:, y, x] = spei_yx['SPEI1'].values
        spei3_tyx_empty[:, y, x] = spei_yx['SPEI3'].values
        spei6_tyx_empty[:, y, x] = spei_yx['SPEI6'].values

#%%      
    #SPI
    
spei_ref:list = [1991, 2020]
distribution:str = 'gamma'
da_spi1 = xr.DataArray(spi1_tyx_empty,
                       dims=['time', 'latitude', 'longitude'],
                       coords=dict(time=ds.coords['time'].values,
                                  
                                   latitude = ds.coords['latitude'].values,
                                   longitude = ds.coords['longitude'].values),
                       
                       name='SPI1')
da_spi1 = da_spi1.astype('float32')  

da_spi3 = xr.DataArray(spi3_tyx_empty,
                       dims=['time', 'latitude', 'longitude'],
                       coords=dict(time=ds.coords['time'].values,
                                   
                                   latitude = ds.coords['latitude'].values,
                                   longitude = ds.coords['longitude'].values),
                       
                       name='SPI3')
da_spi3 = da_spi3.astype('float32') 

da_spi6 = xr.DataArray(spi6_tyx_empty,
                       dims=['time', 'latitude', 'longitude'],
                       coords=dict(time=ds.coords['time'].values,
                                   
                                   latitude = ds.coords['latitude'].values,
                                   longitude = ds.coords['longitude'].values),
                       
                       name='SPI6')
da_spi6 = da_spi6.astype('float32')   
ds_spi_era5land = xr.merge([da_spi1,
                           da_spi3,
                           da_spi6])

ds_spi_era5land.attrs['source_variable_tp'] = 'tp: ERA5 Land Expver1 Monthly Mean Total Precipitation in mm (org tp *1000)'
ds_spi_era5land.attrs['SPI_function'] = f'R package SPEI, distibution={distribution}, reference period {spei_ref[0]}-{spei_ref[1]}'


out_spi = '/home/evelien/Agri-SiSa/data/output/ERA5land_SPI136_1990_2022_NL.nc' 
ds_spi_era5land.to_netcdf(out_spi)

#%%

#SPEI     

spei_ref:list = [1991, 2020]
distribution:str = 'log-logistic'
da_spei1 = xr.DataArray(spei1_tyx_empty,
                       dims=['time', 'latitude', 'longitude'],
                       coords=dict(time=ds.coords['time'].values,
                                   
                                   latitude = ds.coords['latitude'].values,
                                   longitude = ds.coords['longitude'].values),
                       
                       name='SPEI1')
da_spei1 = da_spi1.astype('float32')  

da_spei3 = xr.DataArray(spei3_tyx_empty,
                       dims=['time', 'latitude', 'longitude'],
                       coords=dict(time=ds.coords['time'].values,
                                   
                                   latitude = ds.coords['latitude'].values,
                                   longitude = ds.coords['longitude'].values),
                       
                       name='SPEI3')
da_spei3 = da_spei3.astype('float32') 

da_spei6 = xr.DataArray(spei6_tyx_empty,
                       dims=['time', 'latitude', 'longitude'],
                      coords=dict(time=ds.coords['time'].values,
                                  
                                  latitude = ds.coords['latitude'].values,
                                  longitude = ds.coords['longitude'].values),
                       
                       name='SPEI6')
da_spei6 = da_spei6.astype('float32')   
ds_spei_era5land = xr.merge([da_spei1,
                           da_spei3,
                           da_spei6])


ds_spei_era5land.attrs['source_variable_pev'] = 'pev: ERA5land Monthly Mean Potential Evaporation in mm (org pev * -1000) '
ds_spei_era5land.attrs['SPEI_function'] = f'R package SPEI,  distibution= {distribution}, SPEI reference period {spei_ref[0]}-{spei_ref[1]}' 
ds_spei_era5land.attrs['source_variable_tp'] = 'tp: ERA5 Land Expver1 Monthly Mean Total Precipitation in mm (org tp *1000)'

out_spei = '/home/evelien/Agri-SiSa/data/output/ERA5land_SPEI136_1990_2022_NL.nc' 
ds_spei_era5land.to_netcdf(out_spei)
