from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors

from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma

mask=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/plot/egu/india_mask.nc','r')
ind = mask.variables['MASK'][:,:]
ind= ma.masked_where(ind<=0.0,ind)
area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area']
isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/ric_irr_fert/output/ric_irr_fert.nc','r')
lonisam1=isam.variables['lon'][:]

ind,lona11=shiftgrid(180.5,ind,lonisam1,start=False)
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/rice_AreaYieldProduction.nc','r')
ncvar_m = nclu.variables['riceData'][0,0,:,:]
ncvar_y = nclu.variables['riceData'][0,1,:,:]
ncvar_a = nclu.variables['riceData'][0,4,:,:]
ncvar_p = nclu.variables['riceData'][0,5,:,:]

latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]
lat_new=N.flipud(latnc)

isam=NetCDFFile('../isamhiscru_rice_eguornl_m3.nc','r')
ryield = isam.variables['yield'][:,:,:]
lonisam = isam.variables['lon'][:]
latisam = isam.variables['lat'][:]

lat_new=N.flipud(latnc)
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)

region=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')

latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]
lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonnc,lat_new)

yield_m31 = interp(ncvar_m,lonnc,lat_new,lon,lat,order=1)

yield_m31= ma.masked_where(yield_m31[:,:]<=0.0,yield_m31)

yield_y = interp(ncvar_y,lonnc,lat_new,lon,lat,order=1)
yield_y= ma.masked_where(yield_y[:,:]<=0.0,yield_y)

m3newy=N.zeros((360,720))
m3newp=N.zeros((360,720))
m3newa=N.zeros((360,720))
m3newm=N.zeros((360,720))
for x in range(0,360):
    for y in range(0,720):
        a1=x*6
        a2=(x+1)*6
        b1=y*6
        b2=(y+1)*6
        for j in range(a1,a2):
            for i in range(b1,b2):
                m3newm[x,y]=ncvar_m[j,i]+m3newm[x,y]
                m3newy[x,y]=ncvar_y[j,i]+m3newy[x,y]
                m3newa[x,y]=ncvar_a[j,i]+m3newa[x,y]
                m3newp[x,y]=ncvar_p[j,i]+m3newp[x,y]

        m3newm[x,y]=m3newm[x,y]/36
        m3newy[x,y]=m3newy[x,y]/36

m3newm= ma.masked_where(m3newm[:,:]<=0.0,m3newm)
m3newy= ma.masked_where(m3newy[:,:]<=0.0,m3newy)
m3newa= ma.masked_where(m3newa[:,:]<=0.0,m3newa)
m3newp= ma.masked_where(m3newp[:,:]<=0.0,m3newp)

gg=m3newp/m3newa
gg= ma.masked_where(gg<=0.0,gg)
ryield= ma.masked_where(ryield<=0.0,ryield)

for x in range(0,115):
	for i in range(0,360):
		for j in range(0,720):
			ryield[x,i,j]= ma.masked_where(gg[i,j]<=0.0,ryield[x,i,j])

ryield= ma.masked_where(ind<=0.0,ryield)
gg= ma.masked_where(ind<=0.0,gg)


aa=gg*m3newa/gridarea/0.0001
isam=ryield*m3newa/gridarea/0.0001

ncfile=NetCDFFile('isam_rice_egu_m3.nc','w',format='NETCDF3_64BIT_OFFSET')
ncfile.createDimension('lat', 360)
ncfile.createDimension('lon', 720)
ncfile.createDimension('time', 115)

times = ncfile.createVariable('time', 'f8', ('time',))
latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
maize= ncfile.createVariable('yieldisam', 'f8', ('time','lat','lon'),fill_value=-9999.)
maize1= ncfile.createVariable('yieldm3', 'f8', ('time','lat','lon'),fill_value=-9999.)
area= ncfile.createVariable('gridarea', 'f8', ('lat','lon'),fill_value=-9999.)
m3area= ncfile.createVariable('m3area', 'f8', ('lat','lon'),fill_value=-9999.)

latitudes[:] = latisam
longitudes[:] = lonisam
maize[:]=isam
maize1[:]=aa
area[:]=gridarea
m3area[:]=m3newa


latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'

maize.units = 't/ha'
maize1.units = 't/ha'
area.units = 'm2'
m3area.units = 'ha'

maize.long_name = 'isam production over gridcell area'
maize1.long_name = 'M3 production over gridcell area'

latitudes.long_name = 'latitude'
longitudes.long_name = 'longitude'


ncfile.close()





