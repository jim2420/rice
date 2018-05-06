from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area']

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/ric_irr_fert/output/ric_irr_fert.nc','r')
lonisam1=isam.variables['lon'][:]
ric_i_f=isam.variables['totalyield'][96:103,:,:]
riceb=N.average(ric_i_f,axis=0)

riceb,lona11=shiftgrid(180.5,riceb,lonisam1,start=False)

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


isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/isamhiscru_rice.nc','r')
grow = isam.variables['yield'][96:103,:,:]
lonisam = isam.variables['lon'][:]
latisam = isam.variables['lat'][:]

ryield=N.average(grow,axis=0)

spam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/plot/spamrice_isam.nc','r')
spamy = spam.variables['ricey_total'][:,:]
spamya = spam.variables['riceya_total'][:,:]
spamp = spam.variables['ricep_total'][:,:]
spama = spam.variables['ricea_total'][:,:]
#spamy= ma.masked_where(spamy<=0.0,spamy)
#spamp= ma.masked_where(spamp<=0.0,spamp)
#spama= ma.masked_where(spama<=0.0,spama)
#spamya= ma.masked_where(spamya<=0.0,spamya)

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
        
fig = plt.figure(figsize=(10,20))


ax1 = fig.add_subplot(321)
#ax1.set_title("M3",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')

x,y = map(lon,lat)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
m3newp=maskoceans(x,y,m3newp)
gg=m3newp/m3newa
gg= ma.masked_where(gg<=0.0,gg)
gg= ma.masked_where(spamy<=0.0,gg)
spamy= ma.masked_where(gg<=0.0,spamy)
spamy= ma.masked_where(spamy<=0.0,spamy)
spamy= ma.masked_where(spamy>=10.0**20,spamy)

gg= ma.masked_where(spamy<=0.0,gg)
spamy= ma.masked_where(gg<=0.0,spamy)

cs1 = map.pcolormesh(x,y,gg,cmap=plt.cm.jet,vmin=0,vmax=8)
plt.axis('off')
#cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=18)



ax1 = fig.add_subplot(322)
#map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
spamy=maskoceans(x,y,spamy)

cs1 = map.pcolormesh(x,y,spamy,cmap=plt.cm.jet,vmin=0,vmax=8)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=18)






ax1 = fig.add_subplot(323)
#ax1.set_title("M3",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')

x,y = map(lon,lat)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
m3newp=maskoceans(x,y,m3newp)

cs1 = map.pcolormesh(x,y,gg*m3newa,cmap=plt.cm.jet,vmin=0,vmax=200000)
plt.axis('off')
#cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=18)



ax1 = fig.add_subplot(324)
#ax1.set_title("ISAM",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')


map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
spamp=maskoceans(x,y,spamp)
spamp= ma.masked_where(spamy<=0.0,spamp)

cs1 = map.pcolormesh(x,y,spamp,cmap=plt.cm.jet,vmin=0,vmax=200000)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=18)




ax1 = fig.add_subplot(325)
#ax1.set_title("M3",fontsize=20)
#map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
riceb= ma.masked_where(ryield<=0.0,riceb)

cs1 = map.pcolormesh(x,y,gg*m3newa/gridarea/0.0001,cmap=plt.cm.jet,vmin=0,vmax=2)
plt.axis('off')
#cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=18)



ax1 = fig.add_subplot(326)
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
spamya= ma.masked_where(spamy<=0.0,spamya)

cs1 = map.pcolormesh(x,y,spamya,cmap=plt.cm.jet,vmin=0,vmax=2)
plt.axis('off')
cbar = map.colorbar(cs1,location='right',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=18)



plt.savefig('rice_m3spam.jpg',bbox_inches='tight')

plt.show()


