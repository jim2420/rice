from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma

isam1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/luh2area_850_2015_corrcrop.nc','r')
meareaisam = isam1.variables['fric_tt'][1130:1160,:,:]#1980-2009
meareaisam= ma.masked_where(meareaisam<=0.0,meareaisam)
print meareaisam.shape
fmearea = isam1.variables['fric_tt'][1051,:,:]#1901

fmearea= ma.masked_where(fmearea<=0.0,fmearea)

mask=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
ind = mask.variables['REGION_MASK'][:,:]
ind= ma.masked_where(ind<=0.0,ind)
ind1=N.zeros((30,360,720))
farea=N.zeros((30,360,720))

for i in range(0,30):
	ind1[i,:,:]=ind[:,:]
	farea[i,:,:]=fmearea
meareaisam= ma.masked_where(ind1!=8,meareaisam)
farea= ma.masked_where(ind1!=8,farea)

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area']

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_irr_fert/output/ric_irr_fert_ch4.nc','r')
lonisam1=isam.variables['lon'][:]

tric_i_f=isam.variables['ch4_flux'][79:109,:,:]
tric_i_f= ma.masked_where(ind1!=8,tric_i_f)
tric_i_f=ma.masked_where(tric_i_f<=0.0,tric_i_f)
tric_i_f=ma.filled(tric_i_f, fill_value=0.)



isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_fert/output/ric_fert_ch4.nc','r')
tric_f=isam.variables['ch4_flux'][79:109,:,:]
tric_f= ma.masked_where(ind1!=8,tric_f)
tric_f=ma.masked_where(tric_f<=0.0,tric_f)
tric_f=ma.filled(tric_f, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_cli/output/ric_cli_ch4.nc','r')
tric_fc=isam.variables['ch4_flux'][79:109,:,:]
tric_fc= ma.masked_where(ind1!=8,tric_fc)
tric_fc=ma.masked_where(tric_fc<=0.0,tric_fc)
tric_fc=ma.filled(tric_fc, fill_value=0.)


isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_irr/output/ric_irr_ch4.nc','r')
tric_i=isam.variables['ch4_flux'][79:109,:,:]
tric_i= ma.masked_where(ind1!=8,tric_i)
tric_i=ma.masked_where(tric_i<=0.0,tric_i)
tric_i=ma.filled(tric_i, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_co2/output/ric_co2_ch4.nc','r')
tric_c=isam.variables['ch4_flux'][79:109,:,:]
tric_c= ma.masked_where(ind1!=8,tric_c)
tric_c=ma.masked_where(tric_c<=0.0,tric_c)
tric_c=ma.filled(tric_c, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_ndep/output/ric_ndep_ch4.nc','r')
tric_n=isam.variables['ch4_flux'][79:109,:,:]
tric_n= ma.masked_where(ind1!=8,tric_n)
tric_n=ma.masked_where(tric_n<=0.0,tric_n)
tric_n=ma.filled(tric_n, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium/output/equilibrium_ch4.nc','r')
ric_i_f=isam.variables['ch4_flux'][79:109,:,:]
ric_i_f= ma.masked_where(ind1!=8,ric_i_f)
ric_i_f=ma.masked_where(ric_i_f<=0.0,ric_i_f)
ric_i_f=ma.filled(ric_i_f, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_fert/output/equilibrium_fert_ch4.nc','r')
ric_f=isam.variables['ch4_flux'][79:109,:,:]
ric_f= ma.masked_where(ind1!=8,ric_f)
ric_f=ma.masked_where(ric_f<=0.0,ric_f)
ric_f=ma.filled(ric_f, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_irr/output/equilibrium_irr_ch4.nc','r')
ric_i=isam.variables['ch4_flux'][79:109,:,:]
ric_i= ma.masked_where(ind1!=8,ric_i)
ric_i=ma.masked_where(ric_i<=0.0,ric_i)
ric_i=ma.filled(ric_i, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_co2/output/equilibrium_co2_ch4.nc','r')
ric_c=isam.variables['ch4_flux'][79:109,:,:]
ric_c= ma.masked_where(ind1!=8,ric_c)
ric_c=ma.masked_where(ric_c<=0.0,ric_c)
ric_c=ma.filled(ric_c, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_cli/output/equilibrium_cli_ch4.nc','r')
ric_cli=isam.variables['ch4_flux'][79:109,:,:]
ric_cli= ma.masked_where(ind1!=8,ric_cli)
ric_cli=ma.masked_where(ric_cli<=0.0,ric_cli)
ric_cli=ma.filled(ric_cli, fill_value=0.)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_ndep/output/equilibrium_ndep_ch4.nc','r')
ric_n=isam.variables['ch4_flux'][79:109,:,:]
ric_n= ma.masked_where(ind1!=8,ric_n)

ric_n=ma.masked_where(ric_n<=0.0,ric_n)
ric_n=ma.filled(ric_n, fill_value=0.)



base=N.zeros(30)
base_i=N.zeros(30)#no irrigation

base_c=N.zeros(30)#fixed co2
base_cli=N.zeros(30)#irrigation and fixed climate
base_f=N.zeros(30) #no fertilizer
base_n=N.zeros(30)
tbase=N.zeros(30)
tbase_i=N.zeros(30)#no irrigation
tbase_n=N.zeros(30)
tbase_c=N.zeros(30)#fixed co2
tbase_cli=N.zeros(30)#irrigation and fixed climate
tbase_f=N.zeros(30) #no fertilizer
tbase_lu=N.zeros(30)
base_lu=N.zeros(30)


tbase_lu=N.sum(tric_i_f*16/12*10000*farea,axis=(1,2))
tbase=N.sum(tric_i_f*16/12*10000*meareaisam,axis=(1,2))
tbase_i=N.sum(tric_f*16/12*10000*meareaisam,axis=(1,2))
tbase_f=N.sum(tric_i*16/12*10000*meareaisam,axis=(1,2))
tbase_cli=N.sum(tric_fc*16/12*10000*meareaisam,axis=(1,2))
tbase_c=N.sum(tric_c*16/12*10000*meareaisam,axis=(1,2))
tbase_n=N.sum(tric_n*16/12*10000*meareaisam,axis=(1,2))

base_lu=N.sum(ric_i_f*16/12*10000*farea,axis=(1,2))
base=N.sum(ric_i_f*16/12*10000*meareaisam,axis=(1,2))
base_i=N.sum(ric_i*16/12*10000*meareaisam,axis=(1,2))
base_f=N.sum(ric_f*16/12*10000*meareaisam,axis=(1,2))
base_cli=N.sum(ric_cli*16/12*10000*meareaisam,axis=(1,2))
base_c=N.sum(ric_c*16/12*10000*meareaisam,axis=(1,2))
base_n=N.sum(ric_n*16/12*10000*meareaisam,axis=(1,2))

yirr=(tbase-tbase_i)/(10**12)
yf=(tbase-tbase_f)/(10**12)
ycli=(tbase-tbase_cli)/(10**12)
yc=(tbase-tbase_c)/(10**12)
yn=(tbase-tbase_n)/(10**12)
ylu=(tbase-tbase_lu)/(10**12)

tyirr=(base_i-base)/(10**12)
tyf=(base_f-base)/(10**12)
tycli=(base_cli-base)/(10**12)
tyc=(base_c-base)/(10**12)
tyn=(base_n-base)/(10**12)
tylu=(base-base_lu)/(10**12)

ayirr=(N.average(yirr[0:10]),N.average(yirr[10:20]),N.average(yirr[20:30]))
aycli=(N.average(ycli[0:10]),N.average(ycli[10:20]),N.average(ycli[20:30]))
ayc=(N.average(yc[0:10]),N.average(yc[10:20]),N.average(yc[20:30]))
ayf=(N.average(yf[0:10]),N.average(yf[10:20]),N.average(yf[20:30]))
ayn=(N.average(yn[0:10]),N.average(yn[10:20]),N.average(yn[20:30]))
aylu=(N.average(ylu[0:10]),N.average(ylu[10:20]),N.average(ylu[20:30]))

sayirr=(N.std(yirr[0:10]),N.std(yirr[10:20]),N.std(yirr[20:30]))
saycli=(N.std(ycli[0:10]),N.std(ycli[10:20]),N.std(ycli[20:30]))
sayf=(N.std(yf[0:10]),N.std(yf[10:20]),N.std(yf[20:30]))
sayc=(N.std(yc[0:10]),N.std(yc[10:20]),N.std(yc[20:30]))
sayn=(N.std(yn[0:10]),N.std(yn[10:20]),N.std(yn[20:30]))
saylu=(N.std(ylu[0:10]),N.std(ylu[10:20]),N.std(ylu[20:30]))



tayirr=(N.average(tyirr[0:10]),N.average(tyirr[10:20]),N.average(tyirr[20:30]))
taycli=(N.average(tycli[0:10]),N.average(tycli[10:20]),N.average(tycli[20:30]))
tayc=(N.average(tyc[0:10]),N.average(tyc[10:20]),N.average(tyc[20:30]))
tayf=(N.average(tyf[0:10]),N.average(tyf[10:20]),N.average(tyf[20:30]))
tayn=(N.average(tyn[0:10]),N.average(tyn[10:20]),N.average(tyn[20:30]))
taylu=(N.average(tylu[0:10]),N.average(tylu[10:20]),N.average(tylu[20:30]))


tsayirr=(N.std(tyirr[0:10]),N.std(tyirr[10:20]),N.std(tyirr[20:30]))
tsaycli=(N.std(tycli[0:10]),N.std(tycli[10:20]),N.std(tycli[20:30]))
tsayf=(N.std(tyf[0:10]),N.std(tyf[10:20]),N.std(tyf[20:30]))
tsayc=(N.std(tyc[0:10]),N.std(tyc[10:20]),N.std(tyc[20:30]))
tsayn=(N.std(tyn[0:10]),N.std(tyn[10:20]),N.std(tyn[20:30]))
tsaylu=(N.std(tylu[0:10]),N.std(tylu[10:20]),N.std(tylu[20:30]))


tsbase=(N.average(tbase[0:10]/(10**12)),N.average(tbase[10:20]/(10**12)),N.average(tbase[20:30]/(10**12)))
sbase=(N.std(tbase[0:10]/(10**12)),N.std(tbase[10:20]/(10**12)),N.std(tbase[20:30]/(10**12)))



fig = plt.figure(figsize=(8,4))
n_groups = 3
ax = fig.add_subplot(111)
#plt.ylim(-10,60)
index = N.arange(n_groups)
bar_width = 0.11
opacity = 0.6
#arects21 = plt.bar(0.05+index, tayc, bar_width, yerr=tsayc,hatch='..',
#                     alpha=opacity,color='black',label='CO$_{2}$')
rects3 = plt.bar(0.05+index+bar_width*1, ayc, bar_width, yerr=sayc,
         alpha=opacity,color='black',label='CO$_{2}$')
#rects2 = plt.bar(0.05+index+bar_width*2, taycli, bar_width, yerr=tsaycli,hatch='..',
#         alpha=opacity,color='blue',
#         label='Climate')
rects22 = plt.bar(0.05+index+bar_width*2, aycli, bar_width, yerr=saycli,
         alpha=opacity,color='blue',
         label='Climate')

#rects02 = plt.bar(0.05+index+bar_width*4, tayirr, bar_width, yerr=tsayirr,hatch='..',
#         alpha=opacity,color='red',
#         label='Irrigation')
rects0 = plt.bar(0.05+index+bar_width*3, ayirr, bar_width, yerr=sayirr,
         alpha=opacity,color='red',
         label='Irrigation')

#rects41 = plt.bar(0.05+index+bar_width*6, tayf, bar_width, yerr=tsayf,hatch='..',
#         alpha=opacity,color='green',
#         label='Nitrogen fertilizer')

rects4 = plt.bar(0.05+index+bar_width*4, ayf, bar_width, yerr=sayf,
         alpha=opacity,color='green',
         label='Nitrogen fertilizer')

#rects412 = plt.bar(0.05+index+bar_width*8, taylu, bar_width, yerr=tsaylu,hatch='..',
#         alpha=opacity,color='m',
#         label='Land conversion')

rects42 = plt.bar(0.05+index+bar_width*5, aylu, bar_width, yerr=saylu,
         alpha=opacity,color='m',
         label='Land conversion')

rects = plt.bar(0.05+index+bar_width*6, tsbase, bar_width, yerr=sbase,
         alpha=opacity,color='gray',
         label='Base')
print tsbase,tayc,taycli,tayirr,tayf


plt.ylabel('Effect on CH4 fluxes (TgCH4/yr)',fontsize=18)
plt.tick_params(axis='both',labelsize=18)
plt.xticks(index + bar_width+0.37, ('1980s','1990s','2000s'))
#leg=plt.legend(loc=2)
#leg.get_frame().set_alpha(0.5)

plt.tight_layout()

plt.savefig('ricep_effect_ch4aogsluiitest.png')
plt.show()

