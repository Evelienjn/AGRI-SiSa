#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 15:12:48 2023

@author: evelien

This file contains functions used in the AGRI-SISA project to compute:
    Drought indices
    Reprojections of coordinates
    Data downloads
"""

import rpy2
import xarray as xr
import os
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri


from pyproj import Proj, transform

#%%
# import functions (and packages) from R
r_time_series = robjects.r('ts')
SPEI_package = importr('SPEI')
r_spei_function = robjects.r['spei']
r_spi_function = robjects.r['spi'] 

def SPEI_1d(da_pr:xr.DataArray,
            da_pet:xr.DataArray,
            spei_ref:list = [1991, 2020],
            distribution:str = 'log-logistic') -> xr.Dataset:
    """
    SPEI = standardized preciptation evatranspiration index
    This function computes the SPEI index using the R SPEI package (https://github.com/sbegueria/SPEI)
    for different SPEI timescales
    In:
    da_pr: precipitation monthy means xarray data_array 1D (time axis) in mm
    da_pet: potential evaporation monthly means xarray data_array 1D (time axis) in mm
    spei_ref: refence time period as list containing [start year, end year] 
    distribution: distribution function used to fit the data to compute the SPEI.
    
    Out:
    xarray dataset with dimension time, and containing multiple timeseries:
    depending on how many spei_scales are provided:
      spei_scale = [1,3,6,12,48] returns SPEI1, SPEI3, SPEI6, SPEI12, SPEI48
    
    """
    
    spei_scale = [1,3,6,12,48]

    #create r array precipiation minus potential evaporation
    r_ts_wb = r_time_series(robjects.FloatVector(da_pr-da_pet), start = robjects.IntVector([da_pr.time.dt.year[0].values, da_pr.time.dt.month[0].values]), frequency = 12)
    #initiate empty list to append results to
    ds_spei_list=[]
    for timescale in spei_scale:
        r_spei = r_spei_function(r_ts_wb, scale=timescale, na_rm=True, ref_start=robjects.IntVector([spei_ref[0], 1]), ref_end=robjects.IntVector([spei_ref[1], 12]), distibution=distribution)
        r_spei_values = pandas2ri.ri2py_vector(r_spei.rx2('fitted'))
        da_spei = r_spei_values
        # create xarray data array
        da_spei = xr.DataArray(da_spei,
                           dims=['time'],
                           coords=dict(time=da_pr.coords['time'].values,
                                       longitude=da_pr.coords['latitude'].values,
                                       latitude = da_pr.coords['longitude'].values),
                           name=f"SPEI{timescale}")
        da_spei = da_spei.astype('float32')  
        ds_spei_list.append(da_spei)
    ds_spei_era5 = xr.merge(ds_spei_list)
    ds_spei_era5.attrs = da_pr.attrs
    ds_spei_era5.attrs['source_variable_tp'] = 'tp: ERA5 Monthly Mean Total Precipitation in mm (org tp *1000)'
    ds_spei_era5.attrs['source_variable_pev'] = 'pev: ERA5 Monthly Mean Potential Evaporation in mm (org pev * -1000'
    ds_spei_era5.attrs['SPEI_function'] = f'R package SPEI,  distibution= {distribution}, SPEI reference period {spei_ref[0]}-{spei_ref[1]}' 
    
    return ds_spei_era5



def SPI_1d(da_pr:xr.DataArray,
            spei_ref:list = [1991, 2020],
            distribution:str = 'gamma') -> xr.Dataset:
    """
    SPI = standardized precipitation index
    This function computes the SPI index using the R SPEI package (https://github.com/sbegueria/SPEI)
    for different SPEI timescales
    In:
    da_pr: precipitation monthy means xarray data_array 1D (time axis) in mm
    da_pet: potential evaporation monthly means xarray data_array 1D (time axis) in mm
    spei_ref: refence time period as list containing [start year, end year] 
    distribution: distribution function used to fit the data to compute the SPEI.
    
    Out:
    xarray dataset with dimension time, and containing multiple timeseries:
    depending on how many spei_scales are provided:
      spei_scale = [1,3,6,12,48] returns SPEI1, SPEI3, SPEI6, SPEI12, SPEI48
    
    """
    
    spi_scale = [1,3,6,12,48]
    #initiate empty list to append results to
    ds_spi_list = []
    
    for timescale in spi_scale: 
        #create r array precipiation 
        r_tp= r_time_series(robjects.FloatVector(da_pr), start = robjects.IntVector([da_pr.time.dt.year[0].values, da_pr.time.dt.month[0].values]), frequency = 12)
        r_spi=r_spi_function(r_tp, scale=timescale, na_rm=True, ref_start=robjects.IntVector([spei_ref[0], 1]), ref_end=robjects.IntVector([spei_ref[1], 12]), distibution=distribution)
        r_spi_values = pandas2ri.ri2py_vector(r_spi.rx2('fitted'))
        da_spi = xr.DataArray(r_spi_values,
                           dims=['time'],
                           coords=dict(time=da_pr.coords['time'].values),
                           name=f"SPI{timescale}")
        da_spi = da_spi.astype('float32')  
        ds_spi_list.append(da_spi)
    ds_spi_era5 = xr.merge(ds_spi_list)
    ds_spi_era5.attrs = da_pr.attrs
    ds_spi_era5.attrs['source_variable_tp'] = 'tp: ERA5 Monthly Mean Total Precipitation in mm (org tp *1000)'
    ds_spi_era5.attrs['SPI_function'] = f'R package SPEI, distibution={distribution}, reference period {spei_ref[0]}-{spei_ref[1]}'
    
    return ds_spi_era5


def get_value_at_coords_ds(ds, var, x, y, expver = 1):
    return ds[var].sel(latitude=y, longitude=x, method='nearest')
    
def get_value_at_coords_da(da, x, y, expver = 1):
    return da.sel(latitude=y, longitude=x, method='nearest')


def get_latlon_from_polygon_gdf(polygon):
    x_centroid = polygon.centroid.x
    y_centroid = polygon.centroid.y
    
    lon_polygon, lat_polygon = reproject_coords(proj_in=crs_rdnew,
                                                proj_out=crs_wgs84,
                                                x1=x_centroid,
                                                y1=y_centroid)

    return lat_polygon, lon_polygon

def reproject_coords(proj_in:str,
                     proj_out:str,
                     x1:float,
                     y1:float):
    
    inProj = Proj(init=proj_in)
    outProj = Proj(init=proj_out)
    x2,y2 = transform(inProj,outProj,x1,y1)
    #print('in ', proj_in, ' x1: ', x1, ' y1: ',y1)
    #print('out ', proj_out, ' x2: ', x2, ' y2: ',y2)
    return x2, y2

