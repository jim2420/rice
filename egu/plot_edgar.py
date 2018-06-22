from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma
import matplotlib.colors as colors

#mask=NetCDFFile('india_mask.nc','r')
#ind = mask.variables['MASK'][:,:]
#ind= ma.masked_where(ind<=0.0,ind)
region=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
ind = region.variables['REGION_MASK'][:,:]


area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area']

area=NetCDFFile('edgarch4_1970_2012.nc','r')
ch4obs = area.variables['ratey'][25,:,:]*1000#1995


isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/ric_irr_fert/output/ric_irr_fert.nc','r')
lonisam1=isam.variables['lon'][:]
ric_i_f=isam.variables['totalyield'][96:103,:,:]
riceb=N.average(ric_i_f,axis=0)
riceb,lona11=shiftgrid(180.5,riceb,lonisam1,start=False)
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


isam=NetCDFFile('../../isam_rice_egu_m3_a.nc','r')
grow = isam.variables['yieldisam'][96:103,:,:]
lonisam = isam.variables['lon'][:]
latisam = isam.variables['lat'][:]
ryield1 = isam.variables['yieldm3'][:,:]
#mearea=isam.variables['m3area'][:,:]
ryield=N.average(grow,axis=0)

region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/isamhiscru_rice_ch4aogs_luh2rice1.nc','r')
ch4 = region.variables['ch4'][94,:,:]*16/12#year1995
mearea = region.variables['ricearea'][94,:,:]#year1995

#region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/nonflood/ric_irr_fert/output/ric_irr_fert_ch4.nc','r')
#ch4 = region.variables['ch4_flux'][94,:,:]*16/12#year1995
#ch4,lona11=shiftgrid(180.5,ch4,lonisam1,start=False)




mearea= ma.masked_where(mearea<=0.0,mearea)
gridarea=ma.masked_where(mearea<=0.0,gridarea)
ryield= ma.masked_where(ryield<=0.0,ryield)
ryield1= ma.masked_where(ryield1<=0.0,ryield1)
ryield= ma.masked_where(ryield1<=0.0,ryield)
ryield1= ma.masked_where(ryield<=0.0,ryield1)

region=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
ma1 = region.variables['rice'][96:103,:,:]
ma2 = region.variables['rice_irrig'][96:103,:,:]
ma1=N.average(ma1,axis=0)
ma2=N.average(ma2,axis=0)

maitotal=ma1+ma2
maitotal= ma.masked_where(maitotal[:,:]<=0.0,maitotal)
latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]

lat_new=N.flipud(latnc)
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)



lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonnc,lat_new)

#yield_m31 = interp(ncvar_m,lonnc,lat_new,lon,lat,order=1)

#yield_m31= ma.masked_where(yield_m31[:,:]<=0.0,yield_m31)

#yield_y = interp(ncvar_y,lonnc,lat_new,lon,lat,order=1)
#yield_y= ma.masked_where(yield_y[:,:]<=0.0,yield_y)

cmap = plt.cm.terrain_r

bounds1=[-0.1,0.0,0.01,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
norm1 = colors.BoundaryNorm(bounds1, cmap.N)

#bounds2=[-1,0.0,2,4,6,8,10,12,14,16,18,20,222]
#norm2 = colors.BoundaryNorm(bounds2, cmap.N)
cmap2=colors.ListedColormap(['#7A944D', '#F7CD2B','#FA8D2E', '#E45E39','#E33125','#DE3636','#702D70','#232963'])

cmap2.set_over('gray')
cmap2.set_under('green')

bounds=[0,1,5,10,20,40,60,90,130]
norm2 = colors.BoundaryNorm(bounds, cmap2.N)

 
fig = plt.figure(figsize=(12,12))


ax1 = fig.add_subplot(221)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

x,y = map(lon,lat)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

ch4= ma.masked_where(ind!=8.0,ch4)
ch4obs= ma.masked_where(ind!=8.0,ch4obs)

ax1.set_title("ISAM (gCH4/m2/yr)",fontsize=20)
cs1 = map.pcolormesh(x,y,ch4,cmap=cmap,norm=norm1)
plt.axis('off')



ax1 = fig.add_subplot(222)
ax1.set_title("EDGAR (gCH4/m2/yr)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

cs1 = map.pcolormesh(x,y,ch4obs,cmap=cmap,norm=norm1)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16)


ax1 = fig.add_subplot(223)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')



map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

ch4= ma.masked_where(ind!=8.0,ch4)
ch4obs= ma.masked_where(ind!=8.0,ch4obs)
a1=ch4*gridarea/(10**9)
ax1.set_title("ISAM (GgCH4/gridcell/yr)",fontsize=20)
cs1 = map.pcolormesh(x,y,a1,cmap=cmap2,norm=norm2)
plt.axis('off')



ax1 = fig.add_subplot(224)
ax1.set_title("EDGAR (GgCH4/gridcell/yr)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
a2=ch4obs*gridarea/(10**9)
cs1 = map.pcolormesh(x,y,a2,cmap=cmap2,norm=norm2)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16)

plt.savefig('isam_rice_ch4.jpg',bbox_inches='tight')
print " isam GgCH4",N.sum(a1),"EDGAR",N.sum(a2)
plt.show()


