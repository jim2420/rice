from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma
from statsmodels.stats.weightstats import DescrStatsW
import matplotlib.colors as colors
region=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
ind = region.variables['REGION_MASK'][:,:]
lonisam1=region.variables['lon'][:]
ind,lona11=shiftgrid(180.5,ind,lonisam1,start=False)

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maize_AreaYieldProduction.nc','r')
ncvar_m = nclu.variables['maizeData'][0,0,:,:]
ncvar_y = nclu.variables['maizeData'][0,1,:,:]
ncvar_a = nclu.variables['maizeData'][0,4,:,:]
ncvar_p = nclu.variables['maizeData'][0,5,:,:]

nclu1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/soybean_AreaYieldProduction.nc','r')
ncvar_ms = nclu1.variables['soybeanData'][0,0,:,:]
ncvar_ys = nclu1.variables['soybeanData'][0,1,:,:]
ncvar_as = nclu1.variables['soybeanData'][0,4,:,:]
ncvar_ps = nclu1.variables['soybeanData'][0,5,:,:]

latm3 = nclu1.variables['latitude'][:]
lonm3 = nclu1.variables['longitude'][:]



region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/mai_irr_fert/output/mai_irr_fert.nc','r')
ma1 = region.variables['totalyield'][95:105,:,:]#1996-2005

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/soy_irr_fert/output/soy_irr_fert.nc','r')
ma2 = region1.variables['totalyield'][95:105,:,:]

dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/mai_fert/output/mai_fert.nc','r')
ma1i = dat2.variables['totalyield'][95:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/fixedirr/mai_irr/output/mai_irr.nc','r')
ma1f = dat3.variables['totalyield'][95:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/fixedirr/mai_co2/output/mai_co2.nc','r')
ma1c = dat4.variables['totalyield'][95:105,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/fixedirr/mai_cli/output/mai_cli.nc','r')
ma1cli = dat5.variables['totalyield'][95:105,:,:]

dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/soy_fert/output/soy_fert.nc','r')
ma2i = dat2.variables['totalyield'][95:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/fixedirr/soy_irr/output/soy_irr.nc','r')
ma2f = dat3.variables['totalyield'][95:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/fixedirr/soy_co2/output/soy_co2.nc','r')
ma2c = dat4.variables['totalyield'][95:105,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/his_cru/fixedirr/soy_cli/output/soy_cli.nc','r')
ma2cli = dat5.variables['totalyield'][95:105,:,:]


isam1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/luh2area_850_2015_corrcrop.nc','r')
meareaisam = isam1.variables['fmai_tt'][1146:1156,:,:]#1980-2009
meareaisam= ma.masked_where(meareaisam<=0.0,meareaisam)
ma1=ma1*meareaisam
ma1i=ma1i*meareaisam
ma1f=ma1f*meareaisam
ma1c=ma1c*meareaisam
ma1cli=ma1cli*meareaisam

isam1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/luh2area_850_2015_corrcrop.nc','r')
meareaisam1 = isam1.variables['fsoy_tt'][1146:1156,:,:]#1980-2009
meareaisam1= ma.masked_where(meareaisam1<=0.0,meareaisam1)
ma2=ma2*meareaisam1
ma2i=ma2i*meareaisam1
ma2f=ma2f*meareaisam1
ma2c=ma2c*meareaisam1
ma2cli=ma2cli*meareaisam1



ma1=N.average(ma1,axis=0)
ma2=N.average(ma2,axis=0)
ma1i=N.average(ma1i,axis=0)
ma2i=N.average(ma2i,axis=0)
ma1c=N.average(ma1c,axis=0)
ma2c=N.average(ma2c,axis=0)
ma1f=N.average(ma1f,axis=0)
ma2f=N.average(ma2f,axis=0)
ma1cli=N.average(ma1cli,axis=0)
ma2cli=N.average(ma2cli,axis=0)






isam1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/code/m3_mai.nc','r')
m31 = isam1.variables['m3area'][:,:]
isam1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/maisoy_cheyenne/code/m3_soy.nc','r')
m32 = isam1.variables['m3area'][:,:]


ma1,lona11=shiftgrid(180.5,ma1,lonisam1,start=False)
ma2,lona11=shiftgrid(180.5,ma2,lonisam1,start=False)
ma1i,lona11=shiftgrid(180.5,ma1i,lonisam1,start=False)
ma2i,lona11=shiftgrid(180.5,ma2i,lonisam1,start=False)
ma1f,lona11=shiftgrid(180.5,ma1f,lonisam1,start=False)
ma2f,lona11=shiftgrid(180.5,ma2f,lonisam1,start=False)
ma1c,lona11=shiftgrid(180.5,ma1c,lonisam1,start=False)
ma2c,lona11=shiftgrid(180.5,ma2c,lonisam1,start=False)
ma1cli,lona11=shiftgrid(180.5,ma1cli,lonisam1,start=False)
ma2cli,lona11=shiftgrid(180.5,ma2cli,lonisam1,start=False)


ma1= ma.masked_where(m31<=0.0,ma1)
ma2= ma.masked_where(m32<=0.0,ma2)

ma1i= ma.masked_where(m31<=0.0,ma1i)
ma2i= ma.masked_where(m32<=0.0,ma2i)
ma1f= ma.masked_where(m31<=0.0,ma1f)
ma2f= ma.masked_where(m32<=0.0,ma2f)
ma1cli= ma.masked_where(m31<=0.0,ma1cli)
ma2cli= ma.masked_where(m32<=0.0,ma2cli)
ma1c= ma.masked_where(m31<=0.0,ma1c)
ma2c= ma.masked_where(m32<=0.0,ma2c)



latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]

latm3_new=N.flipud(latm3)
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)

ncvar_ms=N.flipud(ncvar_ms)
ncvar_ys=N.flipud(ncvar_ys)
ncvar_as=N.flipud(ncvar_as)
ncvar_ps=N.flipud(ncvar_ps)



lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonm3,latm3_new)




ma1= ma.masked_where(ind!=8.0,ma1)
ma2= ma.masked_where(ind!=8.0,ma2)

ma1i= ma.masked_where(ind!=8.0,ma1i)
ma2i= ma.masked_where(ind!=8.0,ma2i)
ma1f= ma.masked_where(ind!=8.0,ma1f)
ma2f= ma.masked_where(ind!=8.0,ma2f)
ma1cli= ma.masked_where(ind!=8.0,ma1cli)
ma2cli= ma.masked_where(ind!=8.0,ma2cli)

ma1c= ma.masked_where(ind!=8.0,ma1c)
ma2c= ma.masked_where(ind!=8.0,ma2c)






cmap = plt.cm.bwr
bounds=[-10,-8,-6,-4,-2,0,10,20,30,40,50]

norm = colors.BoundaryNorm(bounds, cmap.N)

fig = plt.figure(figsize=(20,15))

ax1 = fig.add_subplot(421)
#map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')
x,y = map(lona11,lat)

aa=((ma1/gridarea*10000)-(ma1i/gridarea*10000))/(ma1/gridarea*10000)*100
aa= ma.masked_where(aa==0.0,aa)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,aa,cmap=cmap,norm=norm)

plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",ticks=bounds,extend='both')
cbar.ax.tick_params(labelsize=16)



ax1 = fig.add_subplot(422)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')
bb=((ma2/gridarea*10000)-(ma2i/gridarea*10000))/(ma2/gridarea*10000)*100
bb= ma.masked_where(bb==0.0,bb)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,bb,cmap=cmap,norm=norm)
#cs1 = map.pcolormesh(x,y,ma2,cmap=plt.cm.gist_earth,vmin=0,vmax=5)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",ticks=bounds,extend='both')
cbar.ax.tick_params(labelsize=16)

ax1 = fig.add_subplot(423)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

aa=((ma1/gridarea*10000)-(ma1f/gridarea*10000))/(ma1/gridarea*10000)*100
aa= ma.masked_where(aa==0.0,aa)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,aa,cmap=cmap,norm=norm)

plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",ticks=bounds,extend='both')
cbar.ax.tick_params(labelsize=16)



ax1 = fig.add_subplot(424)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')
bb=((ma2/gridarea*10000)-(ma2f/gridarea*10000))/(ma2/gridarea*10000)*100
bb= ma.masked_where(bb==0.0,bb)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,bb,cmap=cmap,norm=norm)
#cs1 = map.pcolormesh(x,y,ma2,cmap=plt.cm.gist_earth,vmin=0,vmax=5)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",ticks=bounds,extend='both')
cbar.ax.tick_params(labelsize=16)

ax1 = fig.add_subplot(425)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

aa=((ma1/gridarea*10000)-(ma1c/gridarea*10000))/(ma1/gridarea*10000)*100
aa= ma.masked_where(aa==0.0,aa)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,aa,cmap=cmap,norm=norm)

plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",ticks=bounds,extend='both')
cbar.ax.tick_params(labelsize=16)



ax1 = fig.add_subplot(426)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')
bb=((ma2/gridarea*10000)-(ma2c/gridarea*10000))/(ma2/gridarea*10000)*100
bb= ma.masked_where(bb==0.0,bb)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,bb,cmap=cmap,norm=norm)
#cs1 = map.pcolormesh(x,y,ma2,cmap=plt.cm.gist_earth,vmin=0,vmax=5)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",ticks=bounds,extend='both')
cbar.ax.tick_params(labelsize=16)

ax1 = fig.add_subplot(427)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

aa=((ma1/gridarea*10000)-(ma1cli/gridarea*10000))/(ma1/gridarea*10000)*100
aa= ma.masked_where(aa==0.0,aa)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,aa,cmap=cmap,norm=norm)

plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",ticks=bounds,extend='both')
cbar.ax.tick_params(labelsize=16)



ax1 = fig.add_subplot(428)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')
bb=((ma2/gridarea*10000)-(ma2cli/gridarea*10000))/(ma2/gridarea*10000)*100
bb= ma.masked_where(bb==0.0,bb)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,bb,cmap=cmap,norm=norm)
#cs1 = map.pcolormesh(x,y,ma2,cmap=plt.cm.gist_earth,vmin=0,vmax=5)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",ticks=bounds,extend='both')
cbar.ax.tick_params(labelsize=16)

plt.savefig('isamdiffbase_maisoy.jpg',bbox_inches='tight')

plt.show()


