from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma
import matplotlib.colors as colors
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import matplotlib.colorbar as colorbar
area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area']
#1997-2003
area=NetCDFFile('c3_flood.nc','r')
flood = N.average(area.variables['flood'][1147:1154,:,:],axis=0)
area=NetCDFFile('c3ann.nc','r')
c3ann = N.average(area.variables['c3ann'][1147:1154,:,:],axis=0)
#print flood.shape
longf = area.variables['lon'][:]
flood,lona11 = shiftgrid(180.5,flood,longf,start=False)
c3ann,lona11 = shiftgrid(180.5,c3ann,longf,start=False)

#print lona11
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



grow= N.zeros((7, 360, 720))
grow1= N.zeros((7, 360, 720))
years = range(1997, 2004)
base = NetCDFFile ("../isamhiscru_rice_ch4.nc", mode='r')
   
ryield=N.average(base.variables['ch4'][96:103,:,:],axis=0)

lonisam=base.variables['lon'][:]
latisam=base.variables['lat'][:]




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



fig = plt.figure(figsize=(15,8))


ax1 = fig.add_subplot(221)
ax1.set_title("ISAM CH4 flux (gCH4/m2/yr)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon,lat)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
gg=m3newp/m3newa
gg= ma.masked_where(gg[:,:]<=0.0,gg)

ryield1=maskoceans(x,y,ryield)
#ryield1= ma.masked_where(ryield1[:,:]<=0.0,ryield1)
ryield1= ma.masked_where(gg<=0.0,ryield1)

cs1 = map.pcolormesh(x,y,ryield1*16/12,cmap=plt.cm.jet,vmin=0,vmax=30)
plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=18)
total=ryield1*m3newa*10000
m3newa[N.isnan(m3newa)]=0
print (N.sum(total))/(10**12),'Tg CH4/yr'


ax1 = fig.add_subplot(222)
ax1.set_title("ISAM CLM rice CH4 (GgCH4/yr)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

cmap=colors.ListedColormap(['#7A944D', '#F7CD2B','#FA8D2E', '#E45E39','#E33125','#DE3636','#702D70','#232963'])

cmap.set_over('gray')
cmap.set_under('green')

bounds=[0,1,5,10,20,40,60,90,130]
norm = colors.BoundaryNorm(bounds, cmap.N)
gridarea=maskoceans(x,y,gridarea)
ryield3=maskoceans(x,y,ryield)
ryield3= ma.masked_where(maitotal<=0.0,ryield3)
aal1=ryield3*gridarea*maitotal*16/12/(10**9)
aal1= ma.masked_where(aal1[:,:]<=0,aal1)
cs1 = map.pcolormesh(x,y,aal1,norm=norm,cmap=cmap,vmin=0,vmax=130)
plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",extend='max')
cbar.ax.tick_params(labelsize=18)


ax1 = fig.add_subplot(223)
ax1.set_title("ISAM LUH2 rice CH4 (GgCH4/yr)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

ryield2=maskoceans(x,y,ryield)
flood=maskoceans(x,y,flood)
flood= ma.masked_where(flood[:,:]<=0.0,flood)
c3ann=maskoceans(x,y,c3ann)
c3ann= ma.masked_where(c3ann[:,:]<=0.0,c3ann)
flood=ma.filled(flood, fill_value=0.)
c3ann=ma.filled(c3ann, fill_value=0.)
#ryield= ma.masked_where(flood[:,:]<=0.0,ryield)
ryield2= ma.masked_where(flood*c3ann<=0.0,ryield2)
#ryield= ma.masked_where(ryield[:,:]<=0.0,ryield)
flood=ma.filled(flood, fill_value=0.)
c3ann=ma.filled(c3ann, fill_value=0.)
flood[N.isnan(flood)]=0
c3ann[N.isnan(c3ann)]=0
c3ann= ma.masked_where(c3ann[:,:]>=2,c3ann)
flood= ma.masked_where(flood[:,:]>=2,flood)
c3ann= ma.masked_where(c3ann[:,:]<=0,c3ann)
flood= ma.masked_where(flood[:,:]<=0,flood)
flood=ma.filled(flood, fill_value=0.)
c3ann=ma.filled(c3ann, fill_value=0.)

ryield2=ma.filled(ryield2, fill_value=0.)



cmap=colors.ListedColormap(['#7A944D', '#F7CD2B','#FA8D2E', '#E45E39','#E33125','#DE3636','#702D70','#232963'])

cmap.set_over('gray')
cmap.set_under('green')

bounds=[0,1,5,10,20,40,60,90,130]
norm = colors.BoundaryNorm(bounds, cmap.N)
gridarea=maskoceans(x,y,gridarea)
aal=ryield2*flood*gridarea*c3ann*16/12/(10**9)
aal= ma.masked_where(aal[:,:]<=0,aal)
cs1 = map.pcolormesh(x,y,aal,norm=norm,cmap=cmap,vmin=0,vmax=130)
plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",extend='max')
cbar.ax.tick_params(labelsize=18)
#ryield=ma.filled(ryield, fill_value=0.)

#total1=ryield*flood*gridarea*c3ann
#print (N.sum(total1))/(10**12),'Tg CH4/yr luh2'

#total1=flood*gridarea*c3ann
#print (N.sum(total1)),'m2 luh2'


ax1 = fig.add_subplot(224)
ax1.set_title("ISAM M3 rice CH4 (GgCH4/yr)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-12, urcrnrlat=55,llcrnrlon=60, urcrnrlon=155, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

#cmap=colors.ListedColormap(['olivedrab', 'gold','darkorange', 'orange','red','crimson','darkmagenta','midnightblue'])
#cmap=colors.ListedColormap([(122/255,148/255,77/255), (247/255,205/255,43/255),(250/255,141/255,46/255), (228/255,94/255,57/255),(227/255,49/255,37/255),(222/255,54/255,54/255),(112/255,45/255,112/255),(35/255,41/255,99/255)])
cmap=colors.ListedColormap(['#7A944D', '#F7CD2B','#FA8D2E', '#E45E39','#E33125','#DE3636','#702D70','#232963'])

cmap.set_over('gray')
cmap.set_under('green')

bounds=[0,1,5,10,20,40,60,90,130]
norm = colors.BoundaryNorm(bounds, cmap.N)

'''
cmap = plt.cm.bwr
norm = colors.BoundaryNorm(bounds, cmap.N)
'''
aab=m3newa*10000*ryield1*16/12/(10**9)
aab= ma.masked_where(aab[:,:]<=0,aab)

cs1 = map.pcolormesh(x,y,aab,norm=norm,cmap=cmap,vmin=0,vmax=130)
plt.axis('off')
'''
cb3 = colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                boundaries=[-10] + bounds + [10],
                                extend='both',
                                extendfrac='auto',
                                ticks=bounds,
                                spacing='uniform',
                                orientation='horizontal')
'''                               
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",extend='max')
cbar.ax.tick_params(labelsize=18)

plt.savefig('isam_ricech4_avg_combine.jpg',bbox_inches='tight')



plt.show()


