from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
mask=NetCDFFile('india_mask.nc','r')
ind = mask.variables['MASK'][:,:]
ind= ma.masked_where(ind<=0.0,ind)

#from pyresample import geometry,image, kd_tree
spam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/plot/spamrice_isam.nc','r')
spamy = spam.variables['ricey_total'][:,:]
spamy[N.isnan(spamy)]=0.
spamy[N.isinf(spamy)]=0.
spamy= ma.masked_where(spamy<=0.0,spamy)
spamy= ma.masked_where(spamy>=10.0**20,spamy)


region=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
maskregion = region.variables['REGION_MASK'][:,:]
latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]

region=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
maskregion = region.variables['REGION_MASK'][:,:]
latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]
#maskregion = ma.masked_not_equal(maskregion,1)
#print maskregion[:,:]

nclu=NetCDFFile('ric_growing_season_dates.nc4','r')
#print(nclu)
## rainfed
#ncvar_har = nclu.variables['harvest_day'][0,:,:]
#ncvar_plt = nclu.variables['planting_day'][0,:,:]
#ncvar_len = nclu.variables['growing_season_length'][0,:,:]

##irrigated
ncvar_har = nclu.variables['harvest_day'][1,:,:]
ncvar_plt = nclu.variables['planting_day'][1,:,:]
ncvar_len = nclu.variables['growing_season_length'][1,:,:]

latnc = nclu.variables['lat'][:]
lonnc = nclu.variables['lon'][:]
#lon,lat = N.meshgrid(lonnc,latnc)
lat_new=N.flipud(latnc)
ncvar_har=N.flipud(ncvar_har)
ncvar_plt=N.flipud(ncvar_plt)
ncvar_len=N.flipud(ncvar_len)

maskus,lona11 = shiftgrid(180.5,maskregion,lonmask,start=False)
lon,lat = N.meshgrid(lonnc,lat_new)

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
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)

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
gg= ma.masked_where(gg[:,:]<=0.0,gg)

grow= N.zeros((7, 360, 720))
begin= N.zeros((7, 360, 720))
end= N.zeros((7, 360, 720))

years = range(1997, 2004)
#years = range(2000,2002)
for i, year in enumerate(years):
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/ric_irr_fert/output/ric_irr_fert.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/ric_irr_fert_long4/output/ric_irr_ferttest.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/ric_irr_fert_long4/output/ric_irr_fert.bgp-yearly_crop_{0}.nc".format(year), mode='r')


    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
   
    yield1 = base.variables["growseason"][:]
    grow[i, :, :] = yield1[3,:,:]
    yielda1 = base.variables["begingrowdate"][:]
    begin[i, :, :] = yielda1[3,:,:]
    yielda1a = base.variables["offset"][:]
    end[i, :, :] = yielda1a[3,:,:]

grow1=N.average(grow,axis=0)
begin1=N.average(begin,axis=0)
end1=N.average(end,axis=0)

grow2,lona11 = shiftgrid(180.5,grow1,lona1,start=False)
begin2,lona11 = shiftgrid(180.5,begin1,lona1,start=False)
end2,lona11 = shiftgrid(180.5,end1,lona1,start=False)

ind,lona11 = shiftgrid(180.5,ind,lona1,start=False)

 
#http://matplotlib.org/basemap/users/mapsetup.html


fig = plt.figure(figsize=(25,15))


#map.drawstates()
#map.drawcountries()
#clevs = N.arange(0,15,1)
#maskus_new1= N.zeros((2160, 4320))

#clevs = [0,1,2,3,4,5,6,7,8,9]
#ncvar_maize[0,1,:,:][N.isnan(ncvar_maize[0,1,:,:])] = -9999

ncvar_plt[N.isnan(ncvar_plt)] = -9999
ncvar_har[N.isnan(ncvar_har)] = -9999
ncvar_len[N.isnan(ncvar_len)] = -9999
#maizemask[N.isnan(maizemask)] = -9999
#maizemask= ma.masked_where(maizemask<=0,maizemask)

ncvar_plt= ma.masked_where(ncvar_plt<=0,ncvar_plt)
ncvar_har= ma.masked_where(ncvar_har<=0,ncvar_har)
ncvar_len= ma.masked_where(ncvar_len<=0,ncvar_len)

ncvar_plt= ma.masked_where(gg<=0,ncvar_plt)
ncvar_har= ma.masked_where(gg<=0,ncvar_har)
ncvar_len= ma.masked_where(gg<=0,ncvar_len)

ncvar_plt= ma.masked_where(maskus!=8,ncvar_plt)
ncvar_har= ma.masked_where(maskus!=8,ncvar_har)
ncvar_len= ma.masked_where(maskus!=8,ncvar_len)

ncvar_plt= ma.masked_where(spamy<=0,ncvar_plt)
ncvar_har= ma.masked_where(spamy<=0,ncvar_har)
ncvar_len= ma.masked_where(spamy<=0,ncvar_len)


grow2= ma.masked_where(ncvar_len<=0,grow2)
begin2= ma.masked_where(ncvar_plt<=0,begin2)
end2= ma.masked_where(ncvar_har<=0,end2)

grow2= ma.masked_where(grow2<=0,grow2)
begin2= ma.masked_where(begin2<=0,begin2)
end2= ma.masked_where(end2<=0,end2)

ncvar_len= ma.masked_where(grow2<=0,ncvar_len)
ncvar_plt= ma.masked_where(begin2<=0,ncvar_plt)
ncvar_har= ma.masked_where(end2<=0,ncvar_har)

#ncvar_len= ma.masked_where(ind<=0,ncvar_len)
#ncvar_plt= ma.masked_where(ind<=0,ncvar_plt)
#ncvar_har= ma.masked_where(ind<=0,ncvar_har)

#grow2= ma.masked_where(ind<=0,grow2)
#begin2= ma.masked_where(ind<=0,begin2)
#end2= ma.masked_where(ind<=0,end2)

#ncvar_maize1= ma.masked_where(ncvar_maize2<0.01,ncvar_maize1)

#print ncvar_maize[0,1,:,:]
#ncvar_maize[0,1,:,:] = ma.masked_where(maskus_new>1.1,ncvar_maize[0,1,:,:])
#print maskus_new

ax1 = fig.add_subplot(331)
ax1.set_title("GGCMI planting date (day of year)",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
x,y = map(lon,lat)

cs = map.pcolormesh(x,y,ncvar_plt,cmap=plt.cm.jet_r,vmin=0,vmax=365)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')

ax1 = fig.add_subplot(334)
ax1.set_title("ISAM planting date (day of year)",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,begin2,cmap=plt.cm.jet_r,vmin=0,vmax=365)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')

ax1 = fig.add_subplot(337)
ax1.set_title("ISAM-GGCMI planting date (days)",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=5, urcrnrlat=35,llcrnrlon=66, urcrnrlon=100, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,begin2-ncvar_plt,cmap=plt.cm.bwr,vmin=-60,vmax=60)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')



ax1 = fig.add_subplot(332)
ax1.set_title("GGCMI harvest date (day of year)",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=5, urcrnrlat=35,llcrnrlon=66, urcrnrlon=100, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
x,y = map(lon,lat)

cs = map.pcolormesh(x,y,ncvar_har,cmap=plt.cm.jet_r,vmin=0,vmax=365)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')

ax1 = fig.add_subplot(335)
ax1.set_title("ISAM harvest date (day of year)",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=5, urcrnrlat=35,llcrnrlon=66, urcrnrlon=100, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,end2,cmap=plt.cm.jet_r,vmin=0,vmax=365)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')

ax1 = fig.add_subplot(338)
ax1.set_title("ISAM-GGCMI harvest date (days)",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=5, urcrnrlat=35,llcrnrlon=66, urcrnrlon=100, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,end2-ncvar_har,cmap=plt.cm.bwr,vmin=-60,vmax=60)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')



ax1 = fig.add_subplot(333)
ax1.set_title("GGCMI growing season length (days)",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=5, urcrnrlat=35,llcrnrlon=66, urcrnrlon=100, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
x,y = map(lon,lat)

cs = map.pcolormesh(x,y,ncvar_len,cmap=plt.cm.jet_r,vmin=0,vmax=365)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')

ax1 = fig.add_subplot(336)
ax1.set_title("ISAM growing season length (days)",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=5, urcrnrlat=35,llcrnrlon=66, urcrnrlon=100, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,grow2,cmap=plt.cm.jet_r,vmin=0,vmax=365)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')

ax1 = fig.add_subplot(339)
ax1.set_title("ISAM-GGCMI growing season length (days)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-15, urcrnrlat=40,llcrnrlon=55, urcrnrlon=145, resolution='c')

#map = Basemap(projection ='cyl', llcrnrlat=5, urcrnrlat=35,llcrnrlon=66, urcrnrlon=100, resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,grow2-ncvar_len,cmap=plt.cm.bwr,vmin=-60,vmax=60)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=16) 
plt.axis('off')

plt.savefig('ricephe_avg_irr_ssaall.jpg',dpi=300,bbox_inches='tight')
plt.show()

