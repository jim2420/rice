from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma

isam1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/luh2area_850_2015_corrcrop.nc','r')
meareaisam = isam1.variables['fric_tt'][1130:1160,:,:]#1980-2009
meareaisam= ma.masked_where(meareaisam<=0.0,meareaisam)
print meareaisam.shape

mask=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
ind = mask.variables['REGION_MASK'][:,:]
ind= ma.masked_where(ind<=0.0,ind)
ind1=N.zeros((30,360,720))
for i in range(0,30):
	ind1[i,:,:]=ind[:,:]

meareaisam= ma.masked_where(ind1!=8,meareaisam)

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area']

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_irr_fert/output/ric_irr_fert_ch4.nc','r')
lonisam1=isam.variables['lon'][:]

tric_i_f=isam.variables['ch4_flux'][79:109,:,:]
tric_i_f= ma.masked_where(ind1!=8,tric_i_f)



isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_fert/output/ric_fert_ch4.nc','r')
tric_f=isam.variables['ch4_flux'][79:109,:,:]
tric_f= ma.masked_where(ind1!=8,tric_f)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_cli/output/ric_cli_ch4.nc','r')
tric_fc=isam.variables['ch4_flux'][79:109,:,:]
tric_fc= ma.masked_where(ind1!=8,tric_fc)


isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_irr/output/ric_irr_ch4.nc','r')
tric_i=isam.variables['ch4_flux'][79:109,:,:]
tric_i= ma.masked_where(ind1!=8,tric_i)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/ric_co2/output/ric_co2_ch4.nc','r')
tric_c=isam.variables['ch4_flux'][79:109,:,:]
tric_c= ma.masked_where(ind1!=8,tric_c)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium/output/equilibrium_ch4.nc','r')
ric_i_f=isam.variables['ch4_flux'][79:109,:,:]
ric_i_f= ma.masked_where(ind1!=8,ric_i_f)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_fert/output/equilibrium_fert_ch4.nc','r')
ric_f=isam.variables['ch4_flux'][79:109,:,:]
ric_f= ma.masked_where(ind1!=8,ric_f)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_irr/output/equilibrium_irr_ch4.nc','r')
ric_i=isam.variables['ch4_flux'][79:109,:,:]
ric_i= ma.masked_where(ind1!=8,ric_i)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_co2/output/equilibrium_co2_ch4.nc','r')
ric_c=isam.variables['ch4_flux'][79:109,:,:]
ric_c= ma.masked_where(ind1!=8,ric_c)

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/equilibrium_cli/output/equilibrium_cli_ch4.nc','r')
ric_cli=isam.variables['ch4_flux'][79:109,:,:]
ric_cli= ma.masked_where(ind1!=8,ric_cli)




base=N.zeros(30)
base_i=N.zeros(30)#no irrigation

base_c=N.zeros(30)#fixed co2
base_cli=N.zeros(30)#irrigation and fixed climate
base_f=N.zeros(30) #no fertilizer

tbase=N.zeros(30)
tbase_i=N.zeros(30)#no irrigation

tbase_c=N.zeros(30)#fixed co2
tbase_cli=N.zeros(30)#irrigation and fixed climate
tbase_f=N.zeros(30) #no fertilizer

tbase=N.sum(tric_i_f*16/12*10000*meareaisam,axis=(1,2))
tbase_i=N.sum(tric_f*16/12*10000*meareaisam,axis=(1,2))
tbase_f=N.sum(tric_i*16/12*10000*meareaisam,axis=(1,2))
tbase_cli=N.sum(tric_fc*16/12*10000*meareaisam,axis=(1,2))
tbase_c=N.sum(tric_c*16/12*10000*meareaisam,axis=(1,2))

base=N.sum(ric_i_f*16/12*10000*meareaisam,axis=(1,2))
base_i=N.sum(ric_i*16/12*10000*meareaisam,axis=(1,2))
base_f=N.sum(ric_f*16/12*10000*meareaisam,axis=(1,2))
base_cli=N.sum(ric_cli*16/12*10000*meareaisam,axis=(1,2))
base_c=N.sum(ric_c*16/12*10000*meareaisam,axis=(1,2))


yirr=(tbase-tbase_i)/tbase*100
yf=(tbase-tbase_f)/tbase*100
ycli=(tbase-tbase_cli)/tbase*100
yc=(tbase-tbase_c)/tbase*100


tyirr=(base_i-base)/tbase*100

tyf=(base_f-base)/tbase*100
tycli=(base_cli-base)/tbase*100
tyc=(base_c-base)/tbase*100

ayirr=(N.average(yirr[0:10]),N.average(yirr[10:20]),N.average(yirr[20:30]))
aycli=(N.average(ycli[0:10]),N.average(ycli[10:20]),N.average(ycli[20:30]))
ayc=(N.average(yc[0:10]),N.average(yc[10:20]),N.average(yc[20:30]))
ayf=(N.average(yf[0:10]),N.average(yf[10:20]),N.average(yf[20:30]))

sayirr=(N.std(yirr[0:10]),N.std(yirr[10:20]),N.std(yirr[20:30]))
saycli=(N.std(ycli[0:10]),N.std(ycli[10:20]),N.std(ycli[20:30]))
sayf=(N.std(yf[0:10]),N.std(yf[10:20]),N.std(yf[20:30]))
sayc=(N.std(yc[0:10]),N.std(yc[10:20]),N.std(yc[20:30]))

tayirr=(N.average(tyirr[0:10]),N.average(tyirr[10:20]),N.average(tyirr[20:30]))
taycli=(N.average(tycli[0:10]),N.average(tycli[10:20]),N.average(tycli[20:30]))
tayc=(N.average(tyc[0:10]),N.average(tyc[10:20]),N.average(tyc[20:30]))
tayf=(N.average(tyf[0:10]),N.average(tyf[10:20]),N.average(tyf[20:30]))

tsayirr=(N.std(tyirr[0:10]),N.std(tyirr[10:20]),N.std(tyirr[20:30]))
tsaycli=(N.std(tycli[0:10]),N.std(tycli[10:20]),N.std(tycli[20:30]))
tsayf=(N.std(tyf[0:10]),N.std(tyf[10:20]),N.std(tyf[20:30]))
tsayc=(N.std(tyc[0:10]),N.std(tyc[10:20]),N.std(tyc[20:30]))



fig = plt.figure(figsize=(10,4))
n_groups = 3
ax = fig.add_subplot(111)
#plt.ylim(-10,60)
index = N.arange(n_groups)
bar_width = 0.1
opacity = 0.6
arects21 = plt.bar(0.1+index, tayc, bar_width, yerr=tsayc,hatch='..',
                     alpha=opacity,color='black',label='CO$_{2}$')
rects3 = plt.bar(0.1+index+bar_width*1, ayc, bar_width, yerr=sayc,
         alpha=opacity,color='black',label='CO$_{2}$')
rects2 = plt.bar(0.1+index+bar_width*2, taycli, bar_width, yerr=tsaycli,hatch='..',
         alpha=opacity,color='blue',
         label='Climate')
rects22 = plt.bar(0.1+index+bar_width*3, aycli, bar_width, yerr=saycli,
         alpha=opacity,color='blue',
         label='Climate')

rects02 = plt.bar(0.1+index+bar_width*4, tayirr, bar_width, yerr=tsayirr,hatch='..',
         alpha=opacity,color='red',
         label='Irrigation')
rects0 = plt.bar(0.1+index+bar_width*5, ayirr, bar_width, yerr=sayirr,
         alpha=opacity,color='red',
         label='Irrigation')

rects41 = plt.bar(0.1+index+bar_width*6, tayf, bar_width, yerr=tsayf,hatch='..',
         alpha=opacity,color='green',
         label='Nitrogen fertilizer')

rects4 = plt.bar(0.1+index+bar_width*7, ayf, bar_width, yerr=sayf,
         alpha=opacity,color='green',
         label='Nitrogen fertilizer')


plt.ylabel('Effect on CH4 fluxes (%)',fontsize=18)
plt.tick_params(axis='both',labelsize=18)
plt.xticks(index + bar_width+0.4, ('1980s','1990s','2000s'))
#leg=plt.legend(loc=2)
#leg.get_frame().set_alpha(0.5)

plt.tight_layout()

plt.savefig('ricech4_effect_aogspecent.png')
plt.show()

