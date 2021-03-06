from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma

mask=NetCDFFile('india_mask.nc','r')
ind = mask.variables['MASK'][:,:]
ind= ma.masked_where(ind<=0.0,ind)


area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area']

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium/output/equilibrium.nc','r')
lonisam1=isam.variables['lon'][:]
ind,lona11=shiftgrid(180.5,ind,lonisam1,start=False)

ric_i_f=isam.variables['totalyield'][79:109,:,:]
ric_i_f,lona11=shiftgrid(180.5,ric_i_f,lonisam1,start=False)
ric_i_f=ma.masked_where(ric_i_f<=0.0,ric_i_f)
ric_i_f=ma.filled(ric_i_f, fill_value=0.)



isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium_irr/output/equilibrium_irr.nc','r')
ric_f=isam.variables['totalyield'][79:109,:,:]
ric_f,lona11=shiftgrid(180.5,ric_f,lonisam1,start=False)

ric_f=ma.masked_where(ric_f<=0.0,ric_f)
ric_f=ma.filled(ric_f, fill_value=0.)


isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/ric_irrfixed_fert/output/ric_irrfixed_fert.nc','r')
#isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/ric_fert_constcli/output/ric_fert_constcli.nc','r')
ric_fc=isam.variables['totalyield'][79:109,:,:]
ric_fc,lona11=shiftgrid(180.5,ric_fc,lonisam1,start=False)
ric_fc=ma.masked_where(ric_fc<=0.0,ric_fc)
ric_fc=ma.filled(ric_fc, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium_fert/output/equilibrium_fert.nc','r')
ric_i=isam.variables['totalyield'][79:109,:,:]
ric_i,lona11=shiftgrid(180.5,ric_i,lonisam1,start=False)
ric_i=ma.masked_where(ric_i<=0.0,ric_i)
ric_i=ma.filled(ric_i, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium_co2/output/equilibrium_co2.nc','r')
ric_c=isam.variables['totalyield'][79:109,:,:]
ric_c,lona11=shiftgrid(180.5,ric_c,lonisam1,start=False)
ric_c=ma.masked_where(ric_c<=0.0,ric_c)
ric_c=ma.filled(ric_c, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium_cli/output/equilibrium_cli.nc','r')
#isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu/ric_irrfixed_fert/output/ric_irrfixed_fert.nc','r')
ric_cli=isam.variables['totalyield'][79:109,:,:]
ric_cli,lona11=shiftgrid(180.5,ric_cli,lonisam1,start=False)
ric_cli=ma.masked_where(ric_cli<=0.0,ric_cli)
ric_cli=ma.filled(ric_cli, fill_value=0.)



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
        
m3newa= ma.masked_where(ind<=0.0,m3newa)
m3newa=ma.filled(m3newa, fill_value=0.)
m3newa[N.isnan(m3newa)] = 0
m3newa[N.isinf(m3newa)] = 0

base=N.zeros(30)
base_i=N.zeros(30)#no irrigation
base_i_cli=N.zeros(30) #no irrigation and fixed climate

base_c=N.zeros(30)#fixed co2
base_cli=N.zeros(30)#irrigation and fixed climate
base_f=N.zeros(30) #no fertilizer



for a in range(0,30):
	for i in range(0,360):
		for j in range(0,720):
			base[a]=base[a]+(m3newa[i,j]*ric_i_f[a,i,j])
                        base_i[a]=base_i[a]+(m3newa[i,j]*ric_f[a,i,j])
                        base_i_cli[a]=base_i_cli[a]+(m3newa[i,j]*ric_fc[a,i,j])

                        base_c[a]=base_c[a]+(m3newa[i,j]*ric_c[a,i,j])
                        base_cli[a]=base_cli[a]+(m3newa[i,j]*ric_cli[a,i,j])
                        base_f[a]=base_f[a]+(m3newa[i,j]*ric_i[a,i,j])


yirr=-(base-base_i)/base*100
#yclinoi=(base-base_i-base_cli+base_i_cli)/(10**6) #varclimate_I-varvlimate_noI-constclim_I+conctcli_noI
#yclinoi=(base-base_i_cli)/(10**6) #varclimate_I-conctcli_noI
yclinoi=(base-base_i_cli)/base*100 #varclimate_I-varclimate_fixedI from conclimate

yf=-(base-base_f)/base*100
ycli=-(base-base_cli)/base*100
yc=-(base-base_c)/base*100



#print yirr
ayirr=(N.average(yirr[0:10]),N.average(yirr[10:20]),N.average(yirr[20:30]))
ayclinoi=(N.average(yclinoi[0:10]),N.average(yclinoi[10:20]),N.average(yclinoi[20:30]))
kayclinoi=(N.average(-yclinoi[0:10]),N.average(-yclinoi[10:20]),N.average(-yclinoi[20:30]))

ayf=(N.average(yf[0:10]),N.average(yf[10:20]),N.average(yf[20:30]))
aycli=(N.average(ycli[0:10]),N.average(ycli[10:20]),N.average(ycli[20:30]))
ayc=(N.average(yc[0:10]),N.average(yc[10:20]),N.average(yc[20:30]))



syirr=(N.std(yirr[0:10]),N.std(yirr[10:20]),N.std(yirr[20:30]))
syclinoi=(N.std(yclinoi[0:10]),N.std(yclinoi[10:20]),N.std(yclinoi[20:30]))

ksyclinoi=(N.std(-yclinoi[0:10]),N.std(-yclinoi[10:20]),N.std(-yclinoi[20:30]))



syf=(N.std(yf[0:10]),N.std(yf[10:20]),N.std(yf[20:30]))
sycli=(N.std(ycli[0:10]),N.std(ycli[10:20]),N.std(ycli[20:30]))
syc=(N.std(yc[0:10]),N.std(yc[10:20]),N.std(yc[20:30]))




fig = plt.figure(figsize=(8,5))
n_groups = 3
ax = fig.add_subplot(111)
#plt.ylim(-10,60)
index = N.arange(n_groups)
bar_width = 0.15
opacity = 0.8
rects3 = plt.bar(0.2+index, ayc, bar_width, yerr=syc,
         alpha=opacity,color='blue',
         label='CO2')
rects2 = plt.bar(0.2+index+bar_width*1, aycli, bar_width, yerr=sycli,
         alpha=opacity,color='black',
         label='Cli')
rects0 = plt.bar(0.2+index+bar_width*2, ayirr, bar_width, yerr=syirr,
         alpha=opacity,color='yellow',
         label='Irri')
rects4 = plt.bar(0.2+index+bar_width*3, ayf, bar_width, yerr=syf,
         alpha=opacity,color='pink',
         label='NF')


plt.ylabel('Effect on rice (%)',fontsize=16)
plt.tick_params(axis='both',labelsize=16)
plt.xticks(index + bar_width+0.35, ('1980s','1990s','2000s'))
leg=plt.legend(loc=2)
leg.get_frame().set_alpha(0.5)

plt.tight_layout()

plt.savefig('rice_effect_egu_directp.png')
plt.show()

