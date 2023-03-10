#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 16:40:27 2023

@author: evelien

compute spi and spei at test coordinates Wagening in RD New (Amersfoort)
"""

#%%
proj = '/home/evelien/Agri-SiSa' 
proj_code = '/home/evelien/Agri-SiSa/code'
import os
print("Current working directory: {0}".format(os.getcwd()))

# Change the current working directory
os.chdir(proj_code)

# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))


#%%
import rpy2
import xarray as xr

import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri


from utils import *
#%%
# import functions (and packages) from R
r_time_series = robjects.r('ts')
SPEI_package = importr('SPEI')
r_spei_function = robjects.r['spei']
r_spi_function = robjects.r['spi'] 

#%%
proj_output_data = os.path.join(proj, 'data', 'output')
spi_out = os.path.join(proj_output_data, 'spi_wageningen.nc')
spei_out=os.path.join(proj_output_data, 'spei_wageningen.nc')
#%%
#read and prepare data from era5
#obtained from data_retrieval/era5_allyears_tp_epot.py
data_path_era5 ='/home/evelien/Agri-SiSa/data/ERA5/ERA5_tp_pev_1959_2022_NL.nc'  

tstart = '1959-01'
tfinal='2021-12'

#select expver=1 (niet erat (prilimilarly), maar era5=validated!)
data = xr.open_dataset(data_path_era5).sel(time=slice(tstart, tfinal), expver=1)

# convert to mm
da_pr = data.tp *1000
da_pr.attrs={'units': 'mm', 'long_name': 'Total precipitation'}

#convert to mm, flip sign because of definition flux away from surface
da_pet =data.pev *-1000
da_pet.attrs={'units': 'mm', 'long_name': 'Potential evaporation'}

#coordinate reference systems
crs_wgs84 ='epsg:4326' 
crs_rdnew = 'epsg:28992' #amersfoort coordinaten
#%%
#test coordinate, wageningen
lat = 51.97; lon=5.66667

#take one gridcell from the era5 grid for the preciptation and potential evaporatin data
da_pr = get_value_at_coords_da(da=da_pr,x=lon,y=lat)
da_pet = get_value_at_coords_da(da=da_pet,x=lon,y=lat)
#%%

spi_wageningen = SPI_1d(da_pr=da_pr)
spei_wageningen = SPEI_1d(da_pr=da_pr, da_pet=da_pet)

#%%
spi_wageningen.to_netcdf(spi_out)
spei_wageningen.to_netcdf(spei_out)