from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors
from statsmodels.stats.weightstats import DescrStatsW
from scipy.stats import ttest_ind
from matplotlib.markers import TICKDOWN
import datetime
from matplotlib.dates import DateFormatter
from scipy import stats

area1=NetCDFFile('/project/projectdirs/m1602/datasets4.full/co2_annual_1765_2016.nc','r')
maifert=area1.variables['co2'][135:240]
#print maifert.shape

mask=NetCDFFile('india_mask.nc','r')
ind = mask.variables['MASK'][:,:]
ind= ma.masked_where(ind<=0.0,ind)


area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area']

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/ric_irr_fert/output/ric_irr_fert.nc','r')
lonisam1=isam.variables['lon'][:]
ind,lona11=shiftgrid(180.5,ind,lonisam1,start=False)
tric_i_f=isam.variables['totalyield'][0:105,:,:]
tric_i_f,lona11=shiftgrid(180.5,tric_i_f,lonisam1,start=False)
tric_i_f=ma.masked_where(tric_i_f<=0.0,tric_i_f)
tric_i_f=ma.filled(tric_i_f, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/ric_irr_fert_constc/output/ric_irr_fert_constc.nc','r')
tric_c=isam.variables['totalyield'][0:105,:,:]
tric_c,lona11=shiftgrid(180.5,tric_c,lonisam1,start=False)
tric_c=ma.masked_where(tric_c<=0.0,tric_c)
tric_c=ma.filled(tric_c, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium/output/equilibrium.nc','r')
ric_i_f=isam.variables['totalyield'][0:105,:,:]
ric_i_f,lona11=shiftgrid(180.5,ric_i_f,lonisam1,start=False)
ric_i_f=ma.masked_where(ric_i_f<=0.0,ric_i_f)
ric_i_f=ma.filled(ric_i_f, fill_value=0.)


isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium_co2/output/equilibrium_co2.nc','r')
ric_c=isam.variables['totalyield'][0:105,:,:]
ric_c,lona11=shiftgrid(180.5,ric_c,lonisam1,start=False)
ric_c=ma.masked_where(ric_c<=0.0,ric_c)
ric_c=ma.filled(ric_c, fill_value=0.)



isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/ric_irr_fert/output/ric_irr_fert_et.nc','r')
ettric_i_f=isam.variables['g_ET'][0:105,3,:,:]
ettric_i_f,lona11=shiftgrid(180.5,ettric_i_f,lonisam1,start=False)
ettric_i_f=ma.masked_where(ettric_i_f<=0.0,ettric_i_f)
ettric_i_f=ma.filled(ettric_i_f, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/ric_irr_fert_constc/output/ric_irr_fert_constc_et.nc','r')
ettric_c=isam.variables['g_ET'][0:105,3,:,:]
ettric_c,lona11=shiftgrid(180.5,ettric_c,lonisam1,start=False)
ettric_c=ma.masked_where(ettric_c<=0.0,ettric_c)
ettric_c=ma.filled(ettric_c, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium/output/equilibrium_et.nc','r')
etric_i_f=isam.variables['g_ET'][0:105,3,:,:]
etric_i_f,lona11=shiftgrid(180.5,etric_i_f,lonisam1,start=False)
etric_i_f=ma.masked_where(etric_i_f<=0.0,etric_i_f)
etric_i_f=ma.filled(etric_i_f, fill_value=0.)


isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/egu_new/equilibrium_co2/output/equilibrium_co2_et.nc','r')
etric_c=isam.variables['g_ET'][0:105,3,:,:]
etric_c,lona11=shiftgrid(180.5,etric_c,lonisam1,start=False)
etric_c=ma.masked_where(etric_c<=0.0,etric_c)
etric_c=ma.filled(etric_c, fill_value=0.)


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
ma1 = region.variables['rice'][0:105,:,:]
ma2 = region.variables['rice_irrig'][0:105,:,:]
#ma1=N.average(ma1,axis=0)
#ma2=N.average(ma2,axis=0)

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


tbase_c=N.zeros(105)#fixed co2
tbase=N.zeros(105)
base_c=N.zeros(105)#fixed co2
base=N.zeros(105)

ettbase_c=N.zeros(105)#fixed co2
ettbase=N.zeros(105)
etbase_c=N.zeros(105)#fixed co2
etbase=N.zeros(105)




for a in range(0,105):
        for i in range(0,360):
                for j in range(0,720):
                        base[a]=base[a]+(m3newa[i,j]*ric_i_f[a,i,j])
                        base_c[a]=base_c[a]+(m3newa[i,j]*ric_c[a,i,j])
                        tbase[a]=tbase[a]+(m3newa[i,j]*tric_i_f[a,i,j])
                        tbase_c[a]=tbase_c[a]+(m3newa[i,j]*tric_c[a,i,j])

                        etbase[a]=etbase[a]+(m3newa[i,j]*etric_i_f[a,i,j])
                        etbase_c[a]=etbase_c[a]+(m3newa[i,j]*etric_c[a,i,j])
                        ettbase[a]=ettbase[a]+(m3newa[i,j]*ettric_i_f[a,i,j])
                        ettbase_c[a]=ettbase_c[a]+(m3newa[i,j]*ettric_c[a,i,j])

yci=((tbase-tbase_c)/tbase*100)
ycd=((base_c-base)/base*100)
etyci=((ettbase-ettbase_c)/ettbase*100)
etycd=((etbase_c-etbase)/etbase*100)



fertr=N.zeros(105)
for i in range(0,105):
	fertr[i]=(maifert[i]-296.0378)/296.0378*100



fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(221)


ax.scatter(fertr,yci,c='blue',alpha=1,label='Direct X Interactive')
plt.xlabel('Relative change in CO2 (%)',fontsize=16)
plt.ylabel('Relative change in yield (%)',fontsize=16)
plt.tick_params(axis='both',labelsize=16)

slope, intercept, r_value, p_value, std_err = stats.linregress(fertr,yci)
line = slope*fertr+intercept
ax.plot(fertr,line,linewidth=2.0)
rr=r_value**2
ax.annotate('y = {:05.3f} x'.format(slope), xy=(0.14, 0.70),xycoords='axes fraction',fontsize=18,color='blue')
ax.annotate('R$^2$ = {:04.2f}'.format(rr), xy=(0.11, 0.62),xycoords='axes fraction',fontsize=18,color='blue')

print 'p',p_value,"r-squared:", r_value**2


ax.scatter(fertr,ycd,c='red',alpha=1,label='Direct')
slope, intercept, r_value, p_value, std_err = stats.linregress(fertr,ycd)
line = slope*fertr+intercept
ax.plot(fertr,line,linewidth=2.0,color='red')
rr=r_value**2
ax.annotate('y = {:05.3f} x'.format(slope), xy=(0.54, 0.15),xycoords='axes fraction',fontsize=18,color='red')
ax.annotate('R$^2$ = {:04.2f}'.format(rr), xy=(0.51, 0.07),xycoords='axes fraction',fontsize=18,color='red')

print 'p',p_value,"r-squared:", r_value**2
ax.legend(loc='2')


ax = fig.add_subplot(223)

ax.scatter(fertr,etyci,c='blue',alpha=1)
plt.xlabel('Relative change in CO2 (%)',fontsize=16)
plt.ylabel('Relative change in ET (%)',fontsize=16)
plt.tick_params(axis='both',labelsize=16)

slope, intercept, r_value, p_value, std_err = stats.linregress(fertr,etyci)
line = slope*fertr+intercept
ax.plot(fertr,line,linewidth=2.0)
rr=r_value**2
ax.annotate('y = {:05.3f} x'.format(slope), xy=(0.14, 0.70),xycoords='axes fraction',fontsize=18,color='blue')
ax.annotate('R$^2$ = {:04.2f}'.format(rr), xy=(0.11, 0.62),xycoords='axes fraction',fontsize=18,color='blue')

print 'p',p_value,"r-squared:", r_value**2


ax.scatter(fertr,etycd,c='red',alpha=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(fertr,etycd)
line = slope*fertr+intercept
ax.plot(fertr,line,linewidth=2.0,color='red')
rr=r_value**2
ax.annotate('y = {:05.3f} x'.format(slope), xy=(0.54, 0.15),xycoords='axes fraction',fontsize=18,color='red')
ax.annotate('R$^2$ = {:04.2f}'.format(rr), xy=(0.51, 0.07),xycoords='axes fraction',fontsize=18,color='red')

print 'p',p_value,"r-squared:", r_value**2



plt.savefig('richis_et_1902_2005_co2.png')
plt.show()
