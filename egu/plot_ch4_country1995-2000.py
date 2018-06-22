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
grow = isam.variables['yieldisam'][96:103,:,:]
lonisam = isam.variables['lon'][:]
latisam = isam.variables['lat'][:]
ryield1 = isam.variables['yieldm3'][:,:]
mearea=isam.variables['m3area'][:,:]
gridarea=isam.variables['gridarea'][:,:]
ryield=N.average(grow,axis=0)



mearea= ma.masked_where(mearea<=0.0,mearea)

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


region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/isamhiscru_rice_ch4egu.nc','r')
ch4 = N.average(region.variables['ch4'][94:100,:,:],axis=0)#year1995-2000




lat_new=N.flipud(latnc)
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)

#ch4= ma.masked_where(ch4<=8.0,ch4)
ch4= ma.masked_where(ind!=8.0,ch4)

ryield= ma.masked_where(ind!=8.0,ryield)
ryield1= ma.masked_where(ind!=8.0,ryield1)
gridarea= ma.masked_where(ind!=8.0,gridarea)
mearea= ma.masked_where(ind!=8.0,mearea)
mearea= ma.masked_where(ch4<=0.0,mearea)
ch4= ma.masked_where(mearea<=0.0,ch4)
#mearea= ma.masked_where(ryield<=0.0,mearea)
#ch4= ma.masked_where(ryield<=0.0,ch4)

lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonnc,lat_new)
isamp=ch4*16/12*mearea*10000
m3p=ryield1*gridarea/10000

name=["Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia1"]
a1=[45100,45600,42901,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,50400]
a2=[45100,45600,42929,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,44600]
g1=[251,306,194];g2=[453,222.6,430,526];g3=[6078,4090,5876,4110];g4=[152,42,157];g5=[1663,1548,767];g6=[0,5]
g7=[298,331,150];g8=[1653,2281,2909,2600];g9=[1166,1299,1328];g10=[744,636,532,566];g11=[1089,2110,1748,1765];g12=[1651,1560,1246,1755];g13=[159,55]
g14=[126,252,127]
yanp=[N.average(g1),N.average(g2),N.average(g3),N.average(g4),N.average(g5),N.average(g6),N.average(g7),N.average(g8),N.average(g9),N.average(g10),N.average(g11),N.average(g12),N.average(g13),N.average(g14)]
yanpp=[N.std(g1),N.std(g2),N.std(g3),N.std(g4),N.std(g5),N.std(g6),N.std(g7),N.std(g8),N.std(g9),N.std(g10),N.std(g11),N.std(g12),N.std(g13),N.std(g14)]
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
	mareamask=mearea
	areak=gridarea
	if idx==2:
		isampmask=ma.masked_where(cou<a1[idx],isampmask)
                mareamask=ma.masked_where(cou<a1[idx],mareamask)
                areak=ma.masked_where(cou<a1[idx],areak)

                isampmask=ma.masked_where(cou>a2[idx],isampmask)
                mareamask=ma.masked_where(cou>a2[idx],mareamask)
	        areak=ma.masked_where(cou>a2[idx],areak)
	elif idx==13:
	        for i in range(0,360):
        	        for j in range(0,720):
                        	if cou[i,j]!=a1[idx] and cou[i,j]!=a2[idx]:
					isampmask[i,j]=0
					areak[i,j]=0
                			mareamask[i,j]=0
	else:
#		print idx
		isampmask=ma.masked_where(cou!=a1[idx],isampmask)	
		mareamask=ma.masked_where(cou!=a1[idx],mareamask)
                areak=ma.masked_where(cou!=a1[idx],areak)
				
	allarea1[idx]=N.sum(mareamask)
	allproisam[idx]=(N.sum(isampmask))/(10**9)
 	allgrid[idx]=N.sum(areak)
print allarea1
allyisam=allproisam/allarea1
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
bar_width = 0.4
opacity = 1.0
rects1 = plt.bar(index+0.1, allproisam, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.1+ bar_width*1,yanp , bar_width,yerr=yanpp,ecolor='black', capsize=10,
                 alpha=opacity,
                 color='r',
                 label='Other studies')
#rects0 = plt.bar(index + 0.18+bar_width*2, yanpp, bar_width,
#                 alpha=opacity,
#                 color='black',
#                 label='ALGAS or UNFCCC')
plt.xlabel('Country',fontsize=35)
plt.ylabel('Methane (Gg CH4 / year)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.11, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=1,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric19952000ch4_count.png')
plt.show()




