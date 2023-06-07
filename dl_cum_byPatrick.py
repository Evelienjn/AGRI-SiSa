import cdsapi
import numpy as np

c = cdsapi.Client()

syear = 2021

years = [str(year) for year in np.arange(syear, syear+10)]
print(years)

fname = 'LUKE_%s_hrly_cum.nc' % years[0]

c.retrieve(
    'reanalysis-era5-land',
    {
        'variable': [
            'potential_evaporation', 'total_precipitation',
        ],
        'year': years,
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12'],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30','31'
        ],
        'time': ['00:00'],
        'area': [
            60.81, 23.4, 60.8,
            23.41,
        ],
        'format': 'netcdf',
    },
    fname)
print(fname, "ready")