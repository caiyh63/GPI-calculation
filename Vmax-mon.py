# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 21:05:42 2022
This code is used to calculate TC ptential maximum intensity (PI)

Use temperature, specific humidity, mean sea surface pressure, and SST
Data from ERA5

Revised 2/23/2023:  Use 1x1 resolution replacing on 2.5x2.5
Revised 11/24/2023: modify the type of output variable 'time' to GREGORIAN calendar.

@author: Yuhao Cai @SYSU
"""
import xarray as xr
import numpy as np
import pandas as pd

path_ERA   = '/home/ys17-23/Extension2/caiyh/ERA5'
path_Had   = '/home/ys17-23/Extension/HadISST1'

f          = xr.open_dataset(path_ERA + '/ERA5_plevel_1x1_1975_2021.nc') 
fsst       = xr.open_dataset(path_Had + '/HadISST_sst_lonfilp.nc') #
fs         = xr.open_dataset(path_ERA + '/ERA5_Surface_1x1_1959_2022.nc') 
lat        = f['latitude']
lon        = f['longitude']

# pecific humidity (kg/kg)
q      = f['q'].loc['1979-01-01':'2021-12-01',::-1,::-1,:].loc[:,1000:100,-60:60,:]     
# air temperature (K)
t      = f['t'].loc['1979-01-01':'2021-12-01',::-1,::-1,:].loc[:,1000:100,-60:60,:]     
# sea surface temperature
sst    = fsst['sst'].loc['1979-01-01':'2022-01-01',::-1,:].loc[:,-60:60,:]            
# mean surface pressure (Pa)
slp    = fs['msl'].loc['1979-01-01':'2021-12-01',::-1,:].loc[:,-60:60,:]   

print(sst)


P        = t.level
PP       = np.array(P)
lat_inWP = q.latitude
lon_inWP = q.longitude
time     = q.time

lon_grid25 = np.array(lon_inWP)
lat_grid25 = np.array(lat_inWP) 

sst25 = sst.interp(longitude=lon_grid25,latitude=lat_grid25, method='linear')  # SST data is still correponsed to atmospheric data

# 检查单位
print(f.r.units)
print(f.q.units)
print(f.t.units)
print(sst25.units)
print(slp.units)


#以下进行单位转换
TT   = np.array(t)
TC   = TT-273.5                                            # degC

qq   = np.array(q)*1000                                    # grams/kg
RR   = qq/(1.0-qq*0.001)                                   # calculate the mix ratio
MSL  = np.array(slp[:,:,:])/100.                         # convert units "hPa"
SSTC = np.array(sst25)

# Import function of PI calculation (download from Kerry Emanuel website)
import sys
sys.path.append("/home/ys17-23/Extension2/caiyh/pyfunc/potential_intensity")
import pi

tlen, ylen, xlen = len(q.time), len(lat_inWP), len(lon_inWP)
# Mytime = np.zeros((tlen))
vmax   = np.zeros((tlen,ylen,xlen))
#print(sst25)
print('PI datas are calculating')

# call PI
for k in range(tlen):
    for j in range(ylen):
    	for i in range(xlen):
        	POTI  = pi.pi(SSTC[k,j,i], MSL[k,j,i], PP, TC[k,:,j,i], RR[k,:,j,i]
        		,CKCD=0.9,ascent_flag=0,diss_flag=1,V_reduc=0.8,miss_handle=1)
        	vmax[k,j,i]  = POTI[0]

print('data is writing')
# write TIME (YYYYMM)
# ii = 0
# for yy in range(1979,2022):
#     for mm in range(1,13):
#         Mytime[ii]= yy*100+mm
#         ii = ii+1


save_path   = '/home/ys17-23/Extension2/caiyh/cal_GPI/ENGPI-ERA5/vpot_1979_2021ttt.nc'
import os
if os.path.isfile(save_path):
	print('File' + save_path + " exist")
	os.remove(save_path)



##=================Create netcdf File===========================================================
import netCDF4 as nc
import dateutil.parser 
f              = nc.Dataset(save_path,'w',format='NETCDF4')
## Create Dimension
f.createDimension('longitude',size=xlen)
f.createDimension('latitude',size=ylen)
f.createDimension('time',None)


## Create Variable
LON            = f.createVariable('longitude','f4',dimensions='longitude')   # LON
LAT            = f.createVariable('latitude','f4',dimensions='latitude')     # LAT
TIME           = f.createVariable('time','int32',dimensions='time')          # TIME
TIME.setncattr('units', 'hours since 1900-01-01 00:00:00.0')
print(TIME.units)

# :datetime64 dtype process
dt = []
for k in time:
    dt.append(dateutil.parser.parse(str(k.values)))
newtime = nc.date2num(dt, TIME.units)

VP             = f.createVariable('vpot','f4',dimensions=('time','latitude','longitude' )) # POTENTIAL INTENSITY

TIME[:]        = newtime
LON[:]         = lon_inWP
LAT[:]         = lat_inWP
VP[...]        = vmax

TIME.long_name = 'time'
TIME.calendar  = 'gregorian'

LON.long_name  = 'longitude,west is negative'
LON.units      = 'degrees_east'

LAT.long_name  = 'latitude,south is negative'
LAT.units      = 'degrees_north'

VP.long_name   = 'VMAX calculated following Emanuel'
VP.units       = 'm/s'

####close file
f.close()
