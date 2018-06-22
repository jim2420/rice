from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma
import matplotlib.colors as colors

mask=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/Ctry_halfdeg.nc','r')
cou = mask.variables['MASK_Country'][:,:]
cou= ma.masked_where(cou<=0.0,cou)

region=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
ind = region.variables['REGION_MASK'][:,:]

area=NetCDFFile('edgarch4_1970_2012.nc','r')
ch4obs = area.variables['amount'][25,:,:]*1000*86400*365#1995


#area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
#gridarea = area.variables['cell_area']

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
gridarea=isam.variables['gridarea'][:,:]
ryield=N.average(grow,axis=0)




ryield= ma.masked_where(ryield<=0.0,ryield)
ryield1= ma.masked_where(ryield1<=0.0,ryield1)

ryield= ma.masked_where(ryield1<=0.0,ryield)
ryield1= ma.masked_where(ryield<=0.0,ryield1)
#mearea= ma.masked_where(ryield<=0.0,mearea)


region=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
ma1 = region.variables['rice'][96:103,:,:]
ma2 = region.variables['rice_irrig'][96:103,:,:]
ma1=N.average(ma1,axis=0)
ma2=N.average(ma2,axis=0)

maitotal=ma1+ma2
maitotal= ma.masked_where(maitotal[:,:]<=0.0,maitotal)
latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]


region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/isamhiscru_rice_ch4aogs_luh2dry1.nc','r')
ch4 = region.variables['ch4'][94,:,:]#year1995
mearea = region.variables['ricearea'][94,:,:]#year1995

mearea= ma.masked_where(mearea<=0.0,mearea)

#region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/nonflood/ric_irr_fert/output/ric_irr_fert_ch4.nc','r')
#ch4 = region.variables['ch4_flux'][94,:,:]#year1995
#ch4,lona11=shiftgrid(180.5,ch4,lonisam1,start=False)


lat_new=N.flipud(latnc)
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)

#ch4= ma.masked_where(ch4<=0.0,ch4)
ch4= ma.masked_where(ind!=8.0,ch4)
ch4obs= ma.masked_where(ind!=8.0,ch4obs)

ryield= ma.masked_where(ind!=8.0,ryield)
ryield1= ma.masked_where(ind!=8.0,ryield1)
gridarea= ma.masked_where(ind!=8.0,gridarea)
mearea= ma.masked_where(ind!=8.0,mearea)
gridarea= ma.masked_where(mearea<=0.0,gridarea)

#mearea= ma.masked_where(ch4<=0.0,mearea)
#ch4= ma.masked_where(mearea<=0.0,ch4)

lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonnc,lat_new)
isamp=ch4*16/12*mearea

#isamp=ch4*16/12*gridarea


#isamp=ch4*16/12*mearea*10000


m3p=ryield1*gridarea/10000

name=["Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia1"]
a1=[45100,45600,42901,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,50400]
a2=[45100,45600,42929,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,44600]
yanp=N.array([194,430,5876,157,1548,5,298,2909,1299,532,1748,1246,55,127])#yan et al. 2003
yanpp=[306,526,4110,42,767,0,150,2600,1328,631,2110,1270,159,252]#yan et al. 2003
area1=N.array([1497000,2162000,42832000,749000,9942000,30000,2059000,11439000,6033000,3620000,8849000,6766000,544000,673000])
faoy=[10.472,14,10.556,10.668,9.94,11.508,15.6058,21.1428,15.6498,33.4565,15.9936,17.6499,9.8251,17.8917]
faop=[156.7438,302.652,4517.968,94.9014,989.199,3.4524,300.2556,2418.4741,944.0959,1257.5265,1457.4889,1194.1243,55.0105,120.3732]
faoa=[2161800,2161800,42800000,889589,9951700,30000,1924000,11438760,6032654,3758691,9112951,6765600,559900,672787]
yyyan=yanp*(10**9)/(area1*10000)
num=14
allarea1=N.zeros((num))
allproisam=N.zeros((num))
allprom3=N.zeros((num))
allyisam=N.zeros((num))
allym3=N.zeros((num))
allgrid=N.zeros((num))
m3grid=N.zeros((num))
isamgrid=N.zeros((num))
ch4obsf=N.zeros((num))
for idx in xrange(num):
	#print a1[idx],a2[idx]
	isampmask=N.zeros((360,720))
	m3pmask=N.zeros((360,720))
        mareamask=N.zeros((360,720))
        areak=N.zeros((360,720))
	ch4obs1=N.zeros((360,720))
	isampmask=isamp
	mareamask=mearea
	areak=gridarea
	ch4obs1=ch4obs
	if idx==2:
		isampmask=ma.masked_where(cou<a1[idx],isampmask)
                mareamask=ma.masked_where(cou<a1[idx],mareamask)
                areak=ma.masked_where(cou<a1[idx],areak)
                ch4obs1=ma.masked_where(cou<a1[idx],ch4obs1)

		isampmask=ma.masked_where(cou>a2[idx],isampmask)
		mareamask=ma.masked_where(cou>a2[idx],mareamask)
		areak=ma.masked_where(cou>a2[idx],areak)
		ch4obs1=ma.masked_where(cou>a2[idx],ch4obs1)
	elif idx==13:
	        for i in range(0,360):
        	        for j in range(0,720):
                        	if cou[i,j]!=a1[idx] and cou[i,j]!=a2[idx]:
					isampmask[i,j]=0
					areak[i,j]=0
                			mareamask[i,j]=0
					ch4obs1[i,j]=0
	else:
#		print idx
		isampmask=ma.masked_where(cou!=a1[idx],isampmask)	
		mareamask=ma.masked_where(cou!=a1[idx],mareamask)
                areak=ma.masked_where(cou!=a1[idx],areak)
                ch4obs1=ma.masked_where(cou!=a1[idx],ch4obs1)
        ch4obsf[idx]=(N.sum(ch4obs1))/(10**9)
	allarea1[idx]=N.sum(mareamask)
	allproisam[idx]=(N.sum(isampmask))/(10**9)
 	allgrid[idx]=N.sum(areak)
#print allgrid,allarea1
allyisam=allproisam*(10**9)/allarea1
#allyisam=allproisam*(10**9)/(allarea1*10000)


isamgrid=allproisam/allgrid*10000


#print allym3,allyisam
#print megrid,isamgrid
#print allprom3,allproisam
newm3y=N.zeros((360,720))
newisamy=N.zeros((360,720))
#newarea=N.zeros((360,720))
newm3p=N.zeros((360,720))
newisamp=N.zeros((360,720))
newm3y1=N.zeros((360,720))
newisamy1=N.zeros((360,720))
for c in range(0,num):
	for i in range(0,360):
		for j in range(0,720):
			if c==13:
                                if cou[i,j]==a1[c] or cou[i,j]==a2[c]:
                                        newisamy[i,j]= allyisam[c]
                                        newisamy1[i,j]= isamgrid[c]
                                        newisamp[i,j]= allproisam[c]

			else:
				if cou[i,j]>=a1[c] and cou[i,j]<=a2[c]:
					newisamy[i,j]= allyisam[c]		
                                	newisamy1[i,j]= isamgrid[c]
                                	newisamp[i,j]= allproisam[c]
				
cmap = plt.cm.terrain_r
bounds=[0.0,1,10,100,1000,5000,10000,50000,100000,200000]
#bounds=[-0.1,0.0,0.1,0.2]
norm2 = colors.BoundaryNorm(bounds, cmap.N)

bounds1=[-0.1,0.0,0.01,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
norm1 = colors.BoundaryNorm(bounds1, cmap.N)


 
fig = plt.figure(figsize=(20,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.2
opacity = 1.0
rects1 = plt.bar(index+0.1, allproisam, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.1+ bar_width*1,yanp , bar_width,
                 alpha=opacity,
                 color='r',
                 label='Yan et al. (2003)')
rects0 = plt.bar(index + 0.1+bar_width*2, faop, bar_width,
                 alpha=opacity,
                 color='black',
                 label='FAO')
rects3 = plt.bar(index + 0.1+bar_width*3, ch4obsf, bar_width,
                 alpha=opacity,
                 color='green',
                 label='EDGAR')
plt.xlabel('Country',fontsize=35)
plt.ylabel('Methane (Gg CH4 / year)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.21, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=1,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric1995ch4_countricedry1.png')
plt.show()

fig = plt.figure(figsize=(20,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.25
opacity = 1.0
rects1 = plt.bar(index+0.1, allyisam, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.1+ bar_width*1,yyyan , bar_width,
                 alpha=opacity,
                 color='r',
                 label='Yan et al. (2003)')
rects0 = plt.bar(index + 0.1+bar_width*2, faoy, bar_width,
                 alpha=opacity,
                 color='black',
                 label='FAO')
plt.xlabel('Country',fontsize=35)
plt.ylabel('Methane EF (g CH4 / m2)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.21, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=1,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric1995ch4_couendry1.png')
plt.show()

fig = plt.figure(figsize=(20,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.25
opacity = 1.0
rects1 = plt.bar(index+0.1, allarea1/10000, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.1+ bar_width*1,area1 , bar_width,
                 alpha=opacity,
                 color='r',
                 label='Yan et al. (2003)')
rects0 = plt.bar(index + 0.1+bar_width*2, faoa, bar_width,
                 alpha=opacity,
                 color='black',
                 label='FAO')
plt.xlabel('Country',fontsize=35)
plt.ylabel('Rice area (ha)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.21, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=1,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)

print allarea1/10000
plt.tight_layout()
plt.savefig('ric1995area.png')
plt.show()

fig = plt.figure(figsize=(20,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.25
opacity = 1.0
rects1 = plt.bar(index+0.1, allyisam*faoa*10000/(10**9), bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.1+ bar_width*1,yanp , bar_width,
                 alpha=opacity,
                 color='r',
                 label='Yan et al. (2003)')
rects0 = plt.bar(index + 0.1+bar_width*2, faop, bar_width,
                 alpha=opacity,
                 color='black',
                 label='FAO')
plt.xlabel('Country',fontsize=35)
plt.ylabel('Emissions (Gg CH4 / year)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.21, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=1,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric1995correctmethanedry1.png')
plt.show()


