from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors

def yieldout(year):
        bb=year-1900
	region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
	if bb <=105:
		maitrop = region1.variables['rice'][bb-1,:,:]
		maitropi=region1.variables['rice_irrig'][bb-1,:,:]
	else:
                maitrop = region1.variables['rice'][104,:,:]
                maitropi=region1.variables['rice_irrig'][104,:,:]
	maitrop=ma.masked_where(maitrop<=0,maitrop)
        maitrop=ma.filled(maitrop, fill_value=0.)
        gridarea = region1.variables['area'][:,:]
	maitropi=ma.masked_where(maitropi<=0,maitropi)
	maitropi=ma.filled(maitropi, fill_value=0.)

	maizetor=maitrop
	maizetoi=maitropi

	maizetrop=maitrop+maitropi

	maizeto = maitrop+maitropi


        isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu/ric_fert/output/ric_fert_ch4.nc','r')
        clmtropf = isam.variables['ch4_flux'][bb-1,:,:]

	clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu/ric_irr_fert/output/ric_irr_fert_ch4.nc','r')
	clmtropfi = clm2.variables['ch4_flux'][bb-1,:,:]

        lonisam=clm2.variables['lon'][:]


        clmtropf,lonisam1 = shiftgrid(180.5,clmtropf,lonisam,start=False)
        clmtropfi,lonisam1 = shiftgrid(180.5,clmtropfi,lonisam,start=False)
        print lonisam1



	clmtropf= ma.masked_where(maitrop<=0,clmtropf)
	clmtropf=ma.filled(clmtropf, fill_value=0.)

	clmtropfi= ma.masked_where(maitropi<=0,clmtropfi)
	clmtropfi=ma.filled(clmtropfi, fill_value=0.)

        yield_clmtf=clmtropf
        yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
        yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)

        yield_clmtfi=clmtropfi
        yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
        yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)


	area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
	gridarea = area.variables['cell_area'][:,:]
	gridlon = area.variables['lon'][:]
	gridlat=area.variables['lat'][:]
	gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)


	lon2,lat2 = N.meshgrid(gridlon,gridlat)


	map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
	x,y = map(lon2,lat2)

	yield_clmtf=maskoceans(x,y,yield_clmtf)
	yield_clmtf = ma.masked_where(maizeto<=0,yield_clmtf)

	yield_clmtfi=maskoceans(x,y,yield_clmtfi)
	yield_clmtfi = ma.masked_where(maizeto<=0,yield_clmtfi)

	clmy=((yield_clmtf*maizetor*gridarea)+(yield_clmtfi*maizetoi*gridarea))/((maizetoi*gridarea)+(maizetor*gridarea))




       
        return clmy

yief= N.zeros((115, 360, 720))
years = range(1901, 2016)
for i, yeara in enumerate(years):


	yie=yieldout(yeara)
        yief[i, :, :] = yie

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)
#print gridlon


ncfile=NetCDFFile('isamhiscru_rice_ch4egu.nc','w',format='NETCDF3_64BIT_OFFSET')
ncfile.createDimension('lat', 360)
ncfile.createDimension('lon', 720)
ncfile.createDimension('time', 115)

times = ncfile.createVariable('time', 'f8', ('time',))
latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
maize= ncfile.createVariable('ch4', 'f8', ('time','lat','lon'),fill_value=-9999.)
latitudes[:] = gridlat
longitudes[:] = gridlon
times[:] = years
maize[:]=yief



latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'

latitudes.long_name = 'latitude'
longitudes.long_name = 'longitude'


ncfile.close()





