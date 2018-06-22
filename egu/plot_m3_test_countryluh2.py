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
#grow = isam.variables['yieldisam'][96:103,:,:]
lonisam = isam.variables['lon'][:]
latisam = isam.variables['lat'][:]
ryield1 = isam.variables['yieldm3'][:,:]
mearea=isam.variables['m3area'][:,:]
gridarea=isam.variables['gridarea'][:,:]
#ryield=N.average(grow,axis=0)

isam1=NetCDFFile('../../isamhiscru_rice_aogsluh2phe.nc','r')
grow = isam1.variables['yield'][96:103,:,:]
ryield=N.average(grow,axis=0)



mearea= ma.masked_where(mearea<=0.0,mearea)

ryield= ma.masked_where(ryield<=0.0,ryield)
ryield1= ma.masked_where(ryield1<=0.0,ryield1)

ryield= ma.masked_where(ryield1<=0.0,ryield)
ryield1= ma.masked_where(ryield<=0.0,ryield1)
mearea= ma.masked_where(ryield<=0.0,mearea)


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



ryield= ma.masked_where(ind!=8.0,ryield)
ryield1= ma.masked_where(ind!=8.0,ryield1)
gridarea= ma.masked_where(ind!=8.0,gridarea)
mearea= ma.masked_where(ind!=8.0,mearea)

lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonnc,lat_new)
isamp=ryield*mearea
m3p=ryield1*gridarea/10000

name=["Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia1"]
a1=[45100,45600,42901,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,50400]
a2=[45100,45600,42929,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,44600]
faop=[4216465,7203900,127464896,2859900,37627500,44300,4026092,51898000,20986900,12389412,25843878,32529500,2201700,2140800]
faoy=[2.7028,3.0312,2.8508,3.4374,3.4836,1.7038,2.1155,4.4007,3.3301,3.0681,2.6128,4.2432,3.0606,3.0640]
faoa=[1560044,2376600,44712000,832000,10801214,26000,1903159,11793000,6302175,4038085,9891200,7666300,719370,698700]
num=14
allarea1=N.zeros((num))
allproisam=N.zeros((num))
allprom3=N.zeros((num))
allyisam=N.zeros((num))
allym3=N.zeros((num))
allgrid=N.zeros((num))
m3grid=N.zeros((num))
isamgrid=N.zeros((num))
for idx in xrange(num):
	#print a1[idx],a2[idx]
	isampmask=N.zeros((360,720))
	m3pmask=N.zeros((360,720))
        mareamask=N.zeros((360,720))
        areak=N.zeros((360,720))
	isampmask=isamp
	m3pmask=m3p
	mareamask=mearea
	areak=gridarea
	if idx==2:
		isampmask=ma.masked_where(cou<a1[idx],isampmask)
                m3pmask=ma.masked_where(cou<a1[idx],m3pmask)
                mareamask=ma.masked_where(cou<a1[idx],mareamask)
                areak=ma.masked_where(cou<a1[idx],areak)

                isampmask=ma.masked_where(cou>a2[idx],isampmask)
                m3pmask=ma.masked_where(cou>a2[idx],m3pmask)
                mareamask=ma.masked_where(cou>a2[idx],mareamask)
	        areak=ma.masked_where(cou>a2[idx],areak)
	elif idx==13:
	        for i in range(0,360):
        	        for j in range(0,720):
                        	if cou[i,j]!=a1[idx] and cou[i,j]!=a2[idx]:
					isampmask[i,j]=0
                			m3pmask[i,j]=0
					areak[i,j]=0
                			mareamask[i,j]=0
	else:
#		print idx
		isampmask=ma.masked_where(cou!=a1[idx],isampmask)	
        	m3pmask=ma.masked_where(cou!=a1[idx],m3pmask)
		mareamask=ma.masked_where(cou!=a1[idx],mareamask)
                areak=ma.masked_where(cou!=a1[idx],areak)
				
	allarea1[idx]=N.sum(mareamask)
	allproisam[idx]=N.sum(isampmask)
	allprom3[idx]=N.sum(m3pmask)
 	allgrid[idx]=N.sum(areak)
#print allgrid,allarea1
allym3=allprom3/allarea1
allyisam=allproisam/allarea1
megrid=allprom3/allgrid*10000
isamgrid=allproisam/allgrid*10000
faoyy=faop/allgrid*10000
print faoyy
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
                                        newm3y[i,j]=allym3[c]
                                        newisamy[i,j]= allyisam[c]
                                        newm3y1[i,j]=megrid[c]
                                        newisamy1[i,j]= isamgrid[c]
                                        newm3p[i,j]=allprom3[c]
                                        newisamp[i,j]= allproisam[c]

			else:
				if cou[i,j]>=a1[c] and cou[i,j]<=a2[c]:
					newm3y[i,j]=allym3[c]
					newisamy[i,j]= allyisam[c]		
                        	        newm3y1[i,j]=megrid[c]
                                	newisamy1[i,j]= isamgrid[c]
                                	newm3p[i,j]=allprom3[c]
                                	newisamp[i,j]= allproisam[c]
				
cmap = plt.cm.terrain_r
bounds=[0.0,1,10,100,1000,5000,10000,50000,100000,200000]
#bounds=[-0.1,0.0,0.1,0.2]
norm2 = colors.BoundaryNorm(bounds, cmap.N)

bounds1=[-0.1,0.0,0.01,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
norm1 = colors.BoundaryNorm(bounds1, cmap.N)


 
fig = plt.figure(figsize=(33,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.21
opacity = 0.8
rects1 = plt.bar(index+0.18, allyisam, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.18+ bar_width*1,allym3 , bar_width,
                 alpha=opacity,
                 color='g',
                 label='M3')
rects0 = plt.bar(index + 0.18+bar_width*2, faoy, bar_width,
                 alpha=opacity,
                 color='r',
                 label='FAO')
plt.xlabel('Country',fontsize=35)
plt.ylabel('Yield (t / ha cropland)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.28, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=2,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric2000yield1_countluh2phe.png')
plt.show()



fig = plt.figure(figsize=(33,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.21
opacity = 0.8
rects1 = plt.bar(index+0.18, allproisam, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.18+ bar_width*1,allprom3 , bar_width,
                 alpha=opacity,
                 color='g',
                 label='M3')
rects0 = plt.bar(index + 0.18+bar_width*2, faop, bar_width,
                 alpha=opacity,
                 color='r',
                 label='FAO')
plt.xlabel('Country',fontsize=35)
plt.ylabel('Production (tonnes)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.28, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
#leg=plt.legend(loc=2,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric2000pro_countluh2phe.png')
plt.show()


fig = plt.figure(figsize=(33,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.21
opacity = 0.8
rects1 = plt.bar(index+0.18, isamgrid, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.18+ bar_width*1,megrid , bar_width,
                 alpha=opacity,
                 color='g',
                 label='M3')
rects0 = plt.bar(index + 0.18+bar_width*2, faoyy, bar_width,
                 alpha=opacity,
                 color='r',
                 label='FAO')
plt.xlabel('Country',fontsize=35)
plt.ylabel('(t grains / ha country)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.28, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
#leg=plt.legend(loc=2,ncol=2, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric2000yield_countluh2phe.png')
plt.show()

