from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma
import matplotlib.colors as colors

mask=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/Ctry_halfdeg.nc','r')
cou1 = mask.variables['MASK_Country'][:,:]
cou1= ma.masked_where(cou1<=0.0,cou1)

region=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
ind = region.variables['REGION_MASK'][:,:]



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
gridarea1=isam.variables['gridarea'][:,:]
ryield=N.average(grow,axis=0)



#mearea= ma.masked_where(mearea<=0.0,mearea)

ryield= ma.masked_where(ryield<=0.0,ryield)
ryield1= ma.masked_where(ryield1<=0.0,ryield1)

ryield= ma.masked_where(ryield1<=0.0,ryield)
ryield1= ma.masked_where(ryield<=0.0,ryield1)
#mearea= ma.masked_where(ryield<=0.0,mearea)

area=NetCDFFile('edgarch4_1970_2012.nc','r')
ch4obs = area.variables['amount'][25:31,:,:]*365*86400/(10**6)#1995-2000


region=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
ma1 = region.variables['rice'][96:103,:,:]
ma2 = region.variables['rice_irrig'][96:103,:,:]
ma1=N.average(ma1,axis=0)
ma2=N.average(ma2,axis=0)

maitotal=ma1+ma2
maitotal= ma.masked_where(maitotal[:,:]<=0.0,maitotal)
latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]


region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/isamhiscru_rice_ch4aogs_luh2flood3.nc','r')
#ch4 = N.average(region.variables['ch4'][94:100,:,:],axis=0)#year1995-2000
#mearea = N.average(region.variables['ricearea'][94:100,:,:],axis=0)#year1995-2000
ch4 = region.variables['ch4'][94:100,:,:]#year1995-2000
mearea = region.variables['ricearea'][94:100,:,:]#year1995-2000
ind1=N.zeros((6,360,720))
cou=N.zeros((6,360,720))
gridarea=N.zeros((6,360,720))

for i in range(0,6):
	ind1[i,:,:]=ind[:,:]
	cou[i,:,:]=cou1[:,:]
	gridarea[i,:,:]=gridarea1[:,:]
lat_new=N.flipud(latnc)
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)
mearea= ma.masked_where(mearea<=0.0,mearea)

#ch4= ma.masked_where(ch4<=8.0,ch4)
ch4= ma.masked_where(ind1!=8.0,ch4)
ch4obs= ma.masked_where(ind1!=8.0,ch4obs)

gridarea= ma.masked_where(ind1!=8.0,gridarea)
mearea= ma.masked_where(ind1!=8.0,mearea)
#mearea= ma.masked_where(ch4<=0.0,mearea)
#ch4= ma.masked_where(mearea<=0.0,ch4)
#mearea= ma.masked_where(ryield<=0.0,mearea)
#ch4= ma.masked_where(ryield<=0.0,ch4)

lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonnc,lat_new)
isamp=ch4*16/12*mearea
#isamp=ch4*16/12*gridarea
m3p=ryield1*gridarea/10000

name=["Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia1"]
a1=[45100,45600,42901,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,50400]
a2=[45100,45600,42929,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,44600]
#fao data 5/25/2018
f1=[156.7438,158.256,157.7439,158.5681,162.4197,163.3678];f2=[302.652,315.154,324.422,339.304,352.156,332.724];f3=[4517.968,4581.304,4588.6721,4729.3308,4767.0896,4719.7987];
f4=[94.9014,70.4172,73.6176,81.8236,92.8969,88.7578];f5=[989.199,1013.88,1020.1422,1005.9119,1064.8677,1073.6407];f6=[3.4524,3.4305,3.2868,3.0914,2.9105,2.9921];
f7=[300.2556,293.233,300.9873,306.2741,324.5136,297.0032];f8=[2418.4741,2446.1647,2355.4335,2480.0927,2529.3563,2493.3704];
f9=[944.0959,902.7399,846.374,854.2272,971.9733,986.2753];f10=[1257.5265,1321.9118,1285.4891,1060.5851,1338.2061,1351.0019];
f11=[1457.4889,1482.1589,1585.412,1521.2345,1594.5491,1581.959];f12=[1194.1243,1236.1665,1253.0928,1297.412,1350.8558,1353.0974];
f13=[55.0105,54.4054,59.0776,60.6735,70.5024,70.6785];f14=[120.3732,122.642,123.6273,120.6625,123.8803,125.0094 ]


g1=[251,306,194];g2=[453,222.6,430,526];g3=[6078,4090,5876,4110,4990];g4=[152,42,157];g5=[1663,1548,767];g6=[0,5]
g7=[298,331,150];g8=[1653,2281,2909,2600,2680];g9=[1166,1299,1328,1310];g10=[744,636,532,566];g11=[1089,2110,1748,1765,1540];g12=[1651,1560,1246,1755];g13=[159,55]
g14=[126,252,127]
num=14
area1=[1497000,2162000,42832000,749000,9942000,30000,2059000,11439000,6033000,3620000,8849000,6766000,544000,673000]#yan et al.

faop=[N.average(f1),N.average(f2),N.average(f3),N.average(f4),N.average(f5),N.average(f6),N.average(f7),N.average(f8),N.average(f9),N.average(f10),N.average(f11),N.average(f12),N.average(f13),N.average(f14)]
faop_std=[N.std(f1),N.std(f2),N.std(f3),N.std(f4),N.std(f5),N.std(f6),N.std(f7),N.std(f8),N.std(f9),N.std(f10),N.std(f11),N.std(f12),N.std(f13),N.std(f14)]

allarea1=N.zeros((num))
allproisam=N.zeros((num))
allprom3=N.zeros((num))
allyisam=N.zeros((num))
allym3=N.zeros((num))
allgrid=N.zeros((num))
m3grid=N.zeros((num))
isamgrid=N.zeros((num))
isamgrid_std=N.zeros((num))
allyisam_std=N.zeros((num))
allproisam_std=N.zeros((num))
ch4obs1=N.zeros((num))
ch4obs1_std=N.zeros((num))
ch4obs1all=N.zeros((num,6))
for idx in xrange(num):
	#print a1[idx],a2[idx]
	isampmask=N.zeros((6,360,720))
	m3pmask=N.zeros((6,360,720))
        mareamask=N.zeros((6,360,720))
        areak=N.zeros((6,360,720))
	isampmask=isamp
	mareamask=mearea
	areak=gridarea
	ch4=ch4obs
	if idx==2:
		isampmask=ma.masked_where(cou<a1[idx],isampmask)
                mareamask=ma.masked_where(cou<a1[idx],mareamask)
                areak=ma.masked_where(cou<a1[idx],areak)
                ch4=ma.masked_where(cou<a1[idx],ch4)

                isampmask=ma.masked_where(cou>a2[idx],isampmask)
                mareamask=ma.masked_where(cou>a2[idx],mareamask)
	        areak=ma.masked_where(cou>a2[idx],areak)
                ch4=ma.masked_where(cou>a2[idx],ch4)
	elif idx==13:
		for u in range(0,6):
		        for i in range(0,360):
	     	   	        for j in range(0,720):
	                        	if cou[u,i,j]!=a1[idx] and cou[u,i,j]!=a2[idx]:
						isampmask[u,i,j]=0
						areak[u,i,j]=0
	                			mareamask[u,i,j]=0
                                                ch4[u,i,j]=0

	else:
		isampmask=ma.masked_where(cou!=a1[idx],isampmask)	
		mareamask=ma.masked_where(cou!=a1[idx],mareamask)
                areak=ma.masked_where(cou!=a1[idx],areak)
                ch4=ma.masked_where(cou!=a1[idx],ch4)
		
	allarea1[idx]=N.average(N.sum(mareamask,axis=(1,2)))
	allproisam[idx]=N.average((N.sum(isampmask,axis=(1,2)))/(10**9))
	allproisam_std[idx]=N.std((N.sum(isampmask,axis=(1,2)))/(10**9))
	allgrid[idx]=N.average(N.sum(areak,axis=(1,2)))
	ch4obs1[idx]=N.average(N.sum(ch4,axis=(1,2)))
        ch4obs1_std[idx]=N.std(N.sum(ch4,axis=(1,2)))

	ch4obs1all[idx,:]=(N.sum(ch4,axis=(1,2)))

	allyisam[idx]=N.average(N.sum(isampmask,axis=(1,2))/(N.sum(mareamask,axis=(1,2))))
	isamgrid[idx]=N.average(N.sum(isampmask,axis=(1,2))/(N.sum(areak,axis=(1,2))*10000))
        allyisam_std[idx]=N.std(N.sum(isampmask,axis=(1,2))/(N.sum(mareamask,axis=(1,2))))
        isamgrid_std[idx]=N.std(N.sum(isampmask,axis=(1,2))/(N.sum(areak,axis=(1,2))*10000))
a1=ch4obs1all[0,:].tolist();a2=ch4obs1all[1,:].tolist();a3=ch4obs1all[2,:].tolist();a4=ch4obs1all[3,:].tolist();a5=ch4obs1all[4,:].tolist();a6=ch4obs1all[5,:].tolist();a7=ch4obs1all[6,:].tolist();a8=ch4obs1all[7,:].tolist();a9=ch4obs1all[8,:].tolist();a10=ch4obs1all[9,:].tolist();a11=ch4obs1all[10,:].tolist();a12=ch4obs1all[11,:].tolist();a13=ch4obs1all[12,:].tolist();a14=ch4obs1all[13,:].tolist()
g1=g1+f1+a1;
g2=g2+f2+a2;g3=g3+f3+a3;g4=g4+f4+a4;g5=g5+f5+a5;g6=g6+f6+a6;g7=g7+f7+a7;g8=g8+f8+a8;g9=g9+f9+a9;g10=g10+f10+a10;g11=g11+f11+a11;g12=g12+f12+a12;g13=g13+f13+a13;g14=g14+f14+a14
yanp=[N.average(g1),N.average(g2),N.average(g3),N.average(g4),N.average(g5),N.average(g6),N.average(g7),N.average(g8),N.average(g9),N.average(g10),N.average(g11),N.average(g12),N.average(g13),N.average(g14)]
yanpp=[N.std(g1),N.std(g2),N.std(g3),N.std(g4),N.std(g5),N.std(g6),N.std(g7),N.std(g8),N.std(g9),N.std(g10),N.std(g11),N.std(g12),N.std(g13),N.std(g14)]

#print allym3,allyisam
#print megrid,isamgrid
newm3y=N.zeros((360,720))
newisamy=N.zeros((360,720))
#newarea=N.zeros((360,720))
newm3p=N.zeros((360,720))
newisamp=N.zeros((360,720))
newm3y1=N.zeros((360,720))
newisamy1=N.zeros((360,720))
#for c in range(0,num):
#	for i in range(0,360):
#		for j in range(0,720):
#			if c==13:
#                                if cou[i,j]==a1[c] or cou[i,j]==a2[c]:
#                                        newisamy[i,j]= allyisam[c]
#                                        newisamy1[i,j]= isamgrid[c]
#                                        newisamp[i,j]= allproisam[c]
#
#			else:
#				if cou[i,j]>=a1[c] and cou[i,j]<=a2[c]:
#					newisamy[i,j]= allyisam[c]		
#                                	newisamy1[i,j]= isamgrid[c]
#                                	newisamp[i,j]= allproisam[c]
#				
cmap = plt.cm.terrain_r
bounds=[0.0,1,10,100,1000,5000,10000,50000,100000,200000]
#bounds=[-0.1,0.0,0.1,0.2]
norm2 = colors.BoundaryNorm(bounds, cmap.N)

bounds1=[-0.1,0.0,0.01,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
norm1 = colors.BoundaryNorm(bounds1, cmap.N)
yanp=N.asarray(yanp)/(10**3)
yanpp=N.asarray(yanpp)/(10**3)
allproisam=N.asarray(allproisam)/(10**3)
allproisam_std=N.asarray(allproisam_std)/(10**3)


 
fig = plt.figure(figsize=(20,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.4
opacity = 1.0
rects1 = plt.bar(index+0.1, allproisam, bar_width,yerr=allproisam_std,ecolor='black', capsize=5,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.1+ bar_width*1,yanp , bar_width,yerr=yanpp,ecolor='black', capsize=5,
                 alpha=opacity,
                 color='r',
                 label='Other studies')
#rects0 = plt.bar(index + 0.1+bar_width*2, ch4obs1, bar_width,yerr=ch4obs1_std,ecolor='black', capsize=5,
#                 alpha=opacity,
#                 color='green',
#                 label='EDGAR')
#rects3 = plt.bar(index + 0.1+bar_width*3, faop, bar_width,yerr=faop_std,ecolor='black', capsize=5,
#                 alpha=opacity,
#                 color='gray',
#                 label='FAO')

plt.xlabel('Country',fontsize=35)
plt.ylabel('Methane (Million t CH4 / year)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.11, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=1,ncol=2, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ricc19952000ch4_countluh2flood3.png')
plt.show()



fig = plt.figure(figsize=(20,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.4
opacity = 1.0
rects1 = plt.bar(index+0.1, allarea1/10000, bar_width,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
rects2 = plt.bar(index +0.1+ bar_width*1,area1 , bar_width,
                 alpha=opacity,
                 color='r',
                 label='Yan et al. (2003)')
#rects0 = plt.bar(index + 0.18+bar_width*2, yanpp, bar_width,
#                 alpha=opacity,
#                 color='black',
#                 label='ALGAS or UNFCCC')
plt.xlabel('Country',fontsize=35)
plt.ylabel('Cultivation area (ha)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.11, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=1,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ricarea.png')
plt.show()

