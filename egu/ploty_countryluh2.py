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




#area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
#gridarea = area.variables['cell_area']

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/ric_irr_fert/output/ric_irr_fert.nc','r')
lonisam1=isam.variables['lon'][:]
ric_i_f=isam.variables['totalyield'][96:103,:,:]
riceb=N.average(ric_i_f,axis=0)


spam=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/spamrice_isam.nc','r')
spamy=spam.variables['ricey_total'][:,:]
spampro=spam.variables['ricep_total'][:,:]
spamarea=spam.variables['ricea_total'][:,:]

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
gridarea1=isam.variables['gridarea'][:,:]
#ryield=N.average(grow,axis=0)



isam1=NetCDFFile('../../isamhiscru_rice_aogsluh2flood3.nc','r')
ryield = isam1.variables['yield'][96:105,:,:]#1997-2005
#ryield=N.average(grow,axis=0)
meareaisam = isam1.variables['ricearea'][96:105,:,:]#1997-2005


meareaisam= ma.masked_where(meareaisam<=0.0,meareaisam)

mearea= ma.masked_where(mearea<=0.0,mearea)

ryield= ma.masked_where(ryield<=0.0,ryield)
ryield1= ma.masked_where(ryield1<=0.0,ryield1)



region=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')

latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]

lat_new=N.flipud(latnc)
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)

ind1=N.zeros((9,360,720))
cou=N.zeros((9,360,720))
gridarea=N.zeros((9,360,720))
for i in range(0,9):
        ind1[i,:,:]=ind[:,:]
        cou[i,:,:]=cou1[:,:]
        gridarea[i,:,:]=gridarea1[:,:]



aa1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/iizumi/iizumi.2013JAN29.rice.1982-2006.30min.nc4','r')
iizumiy=N.zeros((9,360,720))
iizumia=N.zeros((9,360,720))
for i in range(1997,2006):
        xx=i-1982
        j=i-1997
        ya = N.flipud(aa1.variables['yield50'][xx,:,:])#t/ha
        yaa = N.flipud(aa1.variables['area'][xx,:,:])##%
        iizumiy[j,:,:]=ya
        iizumia[j,:,:]=yaa
iizumiy= ma.masked_where(iizumiy<=0.0,iizumiy)
iizumia= ma.masked_where(iizumia<=0.0,iizumia)
iizumiy= ma.masked_where(iizumiy>=(10**18),iizumiy)
iizumia= ma.masked_where(iizumia>=(10**18),iizumia)
iizua=iizumia/100*gridarea1/10000
iizup=iizua*iizumiy

ryield= ma.masked_where(ind1!=8.0,ryield)
ryield1= ma.masked_where(ind!=8.0,ryield1)
gridarea= ma.masked_where(ind1!=8.0,gridarea)
mearea= ma.masked_where(ind!=8.0,mearea)
spamy= ma.masked_where(ind!=8.0,spamy)
spampro= ma.masked_where(ind!=8.0,spampro)
spamarea= ma.masked_where(ind!=8.0,spamarea)
meareaisam= ma.masked_where(ind1!=8.0,meareaisam)
iizua= ma.masked_where(ind1!=8.0,iizua)
iizup= ma.masked_where(ind1!=8.0,iizup)

lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonnc,lat_new)
isamp=ryield*meareaisam
m3p=ryield1*gridarea1/10000

name=["Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia1"]
a1=[45100,45600,42901,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,50400]
a2=[45100,45600,42929,46900,40800,41000,41400,50300,45000,50700,47500,48200,44000,44600]
fa1=[1506340,1514210,1550990,1560044,1516980,1544660,1544660,1559436,1541729];fy1=[2.417,2.4434,2.4722,2.7028,2.7454,2.6754,2.6753,2.8573,2.7825]
fa2=[2317300,2423600,2515400,2376600,2114200,2225000,2460600,2519600,2621400];fy2=[2.8048,2.893,3.0744,3.0312,2.7542,3.0192,2.9551,2.9914,3.1742]
fa3=[43469800,44802300,45160000,44712000,44900000,41176100,42592500,41906700,43659800];fy3=[2.8457,2.8805,2.9782,2.8508,3.1158,2.6163,3.1177,2.9756,3.1537]
fa4=[690079,767000,870800,832000,765040,819590,911440,719690,915260];fy4=[3.2451,3.5102,3.281,3.4374,3.5228,3.4889,3.3696,3.6516,3.5465]
fa5=[10263000,10119838,10712955,10801214,10661000,10771000,10725040,10248102,10524067];fy5=[2.7431,2.9358,3.2139,3.4836,3.402,3.4902,3.5768,3.5359,3.7814]
fa6=[28561,26863,25291,26000,23000,18708,19633,18852,25276];fy6=[1.7214,1.7892,1.8471,1.7038,1.7391,1.9987,2.2871,2.8816,2.6896]
fa7=[1928689,1962566,2079442,1903159,1980295,1994645,2242036,2109050,2414500];fy7=[1.7706,1.7884,1.9433,2.1155,2.0699,1.9164,2.1012,1.9773,2.4793]
fa8=[11140594,11730200,11963204,11793000,11500000,11521166,11477357,11922974,11839060];fy8=[4.4322,4.1974,4.2519,4.4007,4.3879,4.4691,4.5426,4.5365,4.5739]
fa9=[5408224,5458405,6210787,6302175,6411996,6376637,6527585,6807628,7383901];fy9=[3.0308,3.0793,3.2405,3.3301,3.3639,3.3655,3.4883,3.5785,3.6899]
fa10=[3842270,3170042,3999839,4038085,4065441,4046318,4006421,4126645,4070421];fy10=[2.9329,2.6986,2.9468,3.0681,3.1866,3.2797,3.3696,3.513,3.5876]
fa11=[9912790,9511520,9969920,9891200,10125424,9653534,10163878,9992868,10224967];fy11=[2.3788,2.418,2.4244,2.6128,2.8739,2.9338,2.9339,2.8895,2.9974]
fa12=[7099700,7350801,7653600,7666300,7492700,7504300,7452200,7445300,7329200];fy12=[3.8768,4.0007,4.1018,4.2432,4.2853,4.5903,4.6387,4.8553,4.8891]
fa13=[601295,617538,717577,719370,746775,738104,756317,770320,736020];fy13=[2.7607,2.7116,2.9304,3.0606,3.1264,3.2739,3.1403,3.2831,3.489]
fa14=[690975,674404,692389,698700,673600,678500,671800,680700,676200];fy14=[3.0676,2.8829,2.9415,3.064,3.1101,3.2385,3.3596,3.326,3.4221]
fp1=N.multiply(fa1,fy1);fp2=N.multiply(fa2,fy2);fp3=N.multiply(fa3,fy3);fp4=N.multiply(fa4,fy4);fp5=N.multiply(fa5,fy5);fp6=N.multiply(fa6,fy6);fp7=N.multiply(fa7,fy7);
fp8=N.multiply(fa8,fy8);fp9=N.multiply(fa9,fy9);fp10=N.multiply(fa10,fy10);fp11=N.multiply(fa11,fy11);fp12=N.multiply(fa12,fy12);fp13=N.multiply(fa13,fy13);fp14=N.multiply(fa14,fy14)


#fao at 2000
#faop=[4216465,7203900,127464896,2859900,37627500,44300,4026092,51898000,20986900,12389412,25843878,32529500,2201700,2140800]
#faoy=[2.7028,3.0312,2.8508,3.4374,3.4836,1.7038,2.1155,4.4007,3.3301,3.0681,2.6128,4.2432,3.0606,3.0640]
#faoa=[1560044,2376600,44712000,832000,10801214,26000,1903159,11793000,6302175,4038085,9891200,7666300,719370,698700]


num=14
spampp=N.zeros((num))
spamaa=N.zeros((num))

allarea1=N.zeros((num))
allproisam=N.zeros((num))
allprom3=N.zeros((num))
allyisam=N.zeros((num))
allym3=N.zeros((num))
allgrid=N.zeros((num))
m3grid=N.zeros((num))
isamgrid=N.zeros((num))
for idx in xrange(num):
	m3pmask=N.zeros((360,720))
        mareamask=N.zeros((360,720))
        areak=N.zeros((360,720))
        spamp1=N.zeros((360,720))
        spama1=N.zeros((360,720))
	m3pmask=m3p
	mareamask=mearea
	areak=gridarea1
	spama1=spamarea
	spamp1=spampro
	if idx==2:
                m3pmask=ma.masked_where(cou1<a1[idx],m3pmask)
                mareamask=ma.masked_where(cou1<a1[idx],mareamask)
                areak=ma.masked_where(cou1<a1[idx],areak)
                spama1=ma.masked_where(cou1<a1[idx],spama1)
                spamp1=ma.masked_where(cou1<a1[idx],spamp1)

                m3pmask=ma.masked_where(cou1>a2[idx],m3pmask)
                mareamask=ma.masked_where(cou1>a2[idx],mareamask)
	        areak=ma.masked_where(cou1>a2[idx],areak)
                spama1=ma.masked_where(cou1>a2[idx],spama1)
                spamp1=ma.masked_where(cou1>a2[idx],spamp1)

	elif idx==13:
	        for i in range(0,360):
        	        for j in range(0,720):
                        	if cou1[i,j]!=a1[idx] and cou1[i,j]!=a2[idx]:
                			m3pmask[i,j]=0
					areak[i,j]=0
                			mareamask[i,j]=0
              				spama1[i,j]=0
                			spamp1[i,j]=0

	else:
#		print idx
        	m3pmask=ma.masked_where(cou1!=a1[idx],m3pmask)
		mareamask=ma.masked_where(cou1!=a1[idx],mareamask)
                areak=ma.masked_where(cou1!=a1[idx],areak)
                spama1=ma.masked_where(cou1!=a1[idx],spama1)
                spamp1=ma.masked_where(cou1!=a1[idx],spamp1)

	allarea1[idx]=N.sum(mareamask)
	allprom3[idx]=N.sum(m3pmask)
 	allgrid[idx]=N.sum(areak)
	spampp[idx]=N.sum(spamp1)
	spamaa[idx]=N.sum(spama1)
allym3=allprom3/allarea1
megrid=allprom3/allgrid*10000
#faoyy=faop/allgrid*10000
spamyy=spampp/spamaa
spamgrid=spampp/allgrid*10000


faoyy1=fp1/allgrid[0]*10000;faoyy2=fp2/allgrid[1]*10000;faoyy3=fp3/allgrid[2]*10000;faoyy4=fp4/allgrid[3]*10000;faoyy5=fp5/allgrid[4]*10000;faoyy6=fp6/allgrid[5]*10000;
faoyy7=fp7/allgrid[6]*10000;faoyy8=fp8/allgrid[7]*10000;faoyy9=fp9/allgrid[8]*10000;faoyy10=fp10/allgrid[9]*10000;faoyy11=fp11/allgrid[10]*10000;faoyy12=fp12/allgrid[11]*10000;faoyy13=fp13/allgrid[12]*10000;faoyy14=fp14/allgrid[13]*10000;



iizuy1=N.zeros((num))
iizuyg1=N.zeros((num))
iizuy1_std=N.zeros((num))
iizuyg1_std=N.zeros((num))

iizua1=N.zeros((num))
iizup1=N.zeros((num))
iizup1_std=N.zeros((num))

iially=N.zeros((num,9))
iiallp=N.zeros((num,9))
iiallg=N.zeros((num,9))

luharea=N.zeros((num))
allyisam_std=N.zeros((num))
allproisam_std=N.zeros((num))
isamgrid_std=N.zeros((num))

for idx in xrange(num):
        #print a1[idx],a2[idx]
        isampmask=N.zeros((9,360,720))
        isamarea=N.zeros((9,360,720))
        areak=N.zeros((9,360,720))
        aiizu=N.zeros((9,360,720))
        piizu=N.zeros((9,360,720))
	aiizu=iizua
	piizu=iizup
        areak=gridarea
        isampmask=isamp
        isamarea=meareaisam
        if idx==2:
                isampmask=ma.masked_where(cou<a1[idx],isampmask)
                isamarea=ma.masked_where(cou<a1[idx],isamarea)
                areak=ma.masked_where(cou<a1[idx],areak)
                aiizu=ma.masked_where(cou<a1[idx],aiizu)
                piizu=ma.masked_where(cou<a1[idx],piizu)

                areak=ma.masked_where(cou>a2[idx],areak)
                isampmask=ma.masked_where(cou>a2[idx],isampmask)
                isamarea=ma.masked_where(cou>a2[idx],isamarea)
                aiizu=ma.masked_where(cou>a2[idx],aiizu)
                piizu=ma.masked_where(cou>a2[idx],piizu)

        elif idx==13:
		for c in range(0,9):
	                for i in range(0,360):
	                        for j in range(0,720):
	                                if cou[c,i,j]!=a1[idx] and cou[c,i,j]!=a2[idx]:
	                                        isampmask[c,i,j]=0
	                                        isamarea[c,i,j]=0
                                                areak[c,i,j]=0
						aiizu[c,i,j]=0
						piizu[c,i,j]=0
        else:
                isampmask=ma.masked_where(cou!=a1[idx],isampmask)
                isamarea=ma.masked_where(cou!=a1[idx],isamarea)
                areak=ma.masked_where(cou!=a1[idx],areak)
                aiizu=ma.masked_where(cou!=a1[idx],aiizu)
                piizu=ma.masked_where(cou!=a1[idx],piizu)
	iizua1[idx]=N.average(N.sum(aiizu,axis=(1,2)))
        luharea[idx]=N.average(N.sum(isamarea,axis=(1,2)))
        iizup1[idx]=N.average((N.sum(piizu,axis=(1,2))))
        iizup1_std[idx]=N.std((N.sum(piizu,axis=(1,2))))

	iizuy1[idx]=N.average(N.sum(piizu,axis=(1,2))/(N.sum(aiizu,axis=(1,2))))
        iizuy1_std[idx]=N.std(N.sum(piizu,axis=(1,2))/(N.sum(aiizu,axis=(1,2))))
        iizuyg1[idx]=N.average((N.sum(piizu,axis=(1,2))/(N.sum(areak,axis=(1,2)))*10000))
        iizuyg1_std[idx]=N.std((N.sum(piizu,axis=(1,2))/(N.sum(areak,axis=(1,2)))*10000))

        allproisam[idx]=N.average((N.sum(isampmask,axis=(1,2))))
        allyisam[idx]=N.average(N.sum(isampmask,axis=(1,2))/(N.sum(isamarea,axis=(1,2))))
        allproisam_std[idx]=N.std((N.sum(isampmask,axis=(1,2))))
        allyisam_std[idx]=N.std(N.sum(isampmask,axis=(1,2))/(N.sum(isamarea,axis=(1,2))))
        isamgrid[idx]=N.average((N.sum(isampmask,axis=(1,2))/(N.sum(areak,axis=(1,2)))*10000))
        isamgrid_std[idx]=N.std((N.sum(isampmask,axis=(1,2))/(N.sum(areak,axis=(1,2)))*10000))

        iiallp[idx,:]=(N.sum(piizu,axis=(1,2)))
        iially[idx,:]=((N.sum(piizu,axis=(1,2)))/(N.sum(aiizu,axis=(1,2))))
        iiallg[idx,:]=((N.sum(piizu,axis=(1,2)))/(N.sum(areak,axis=(1,2)))*10000)


newm3y=N.zeros((360,720))
newisamy=N.zeros((360,720))
#newarea=N.zeros((360,720))
newm3p=N.zeros((360,720))
newisamp=N.zeros((360,720))
newm3y1=N.zeros((360,720))
newisamy1=N.zeros((360,720))


for i in range(1,15):
        kc=([allym3[i-1],spamyy[i-1]])
        locals()['fy{0}'.format(i)]=N.concatenate([locals()['fy{0}'.format(i)],kc])
        locals()['fy{0}'.format(i)]=N.concatenate([locals()['fy{0}'.format(i)],iially[i-1,:]])
    #    print kc,locals()['fy{0}'.format(i)],allym3[i-1],iially[i-1,:]
        locals()['fy{0}'.format(i)]=N.asarray(locals()['fy{0}'.format(i)])

        kc=([megrid[i-1],spamgrid[i-1]])
        locals()['faoyy{0}'.format(i)]=N.concatenate([locals()['faoyy{0}'.format(i)],kc])
        locals()['faoyy{0}'.format(i)]=N.concatenate([locals()['faoyy{0}'.format(i)],iiallg[i-1,:]])
        locals()['faoyy{0}'.format(i)]=N.asarray(locals()['faoyy{0}'.format(i)])

	kc=([allprom3[i-1],spampp[i-1]])
        locals()['fp{0}'.format(i)]=N.concatenate([locals()['fp{0}'.format(i)],kc])
        locals()['fp{0}'.format(i)]=N.concatenate([locals()['fp{0}'.format(i)],iiallp[i-1,:]])
        locals()['fp{0}'.format(i)]=N.asarray(locals()['fp{0}'.format(i)])
        print kc,locals()['fp{0}'.format(i)],allprom3[i-1],iiallp[i-1,:]

print fp10
print fy10
print faoyy10
				
faop=[N.average(fp1),N.average(fp2),N.average(fp3),N.average(fp4),N.average(fp5),N.average(fp6),N.average(fp7),N.average(fp8),N.average(fp9),N.average(fp10),N.average(fp11),N.average(fp12),N.average(fp13),N.average(fp14)]
faop_std=[N.std(fp1),N.std(fp2),N.std(fp3),N.std(fp4),N.std(fp5),N.std(fp6),N.std(fp7),N.std(fp8),N.std(fp9),N.std(fp10),N.std(fp11),N.std(fp12),N.std(fp13),N.std(fp14)]

faoy=[N.average(fy1),N.average(fy2),N.average(fy3),N.average(fy4),N.average(fy5),N.average(fy6),N.average(fy7),N.average(fy8),N.average(fy9),N.average(fy10),N.average(fy11),N.average(fy12),N.average(fy13),N.average(fy14)]
faoy_std=[N.std(fy1),N.std(fy2),N.std(fy3),N.std(fy4),N.std(fy5),N.std(fy6),N.std(fy7),N.std(fy8),N.std(fy9),N.std(fy10),N.std(fy11),N.std(fy12),N.std(fy13),N.std(fy14)]

faoyy=[N.average(faoyy1),N.average(faoyy2),N.average(faoyy3),N.average(faoyy4),N.average(faoyy5),N.average(faoyy6),N.average(faoyy7),N.average(faoyy8),N.average(faoyy9),N.average(faoyy10),N.average(faoyy11),N.average(faoyy12),N.average(faoyy13),N.average(faoyy14)]
faoyy_std=[N.std(faoyy1),N.std(faoyy2),N.std(faoyy3),N.std(faoyy4),N.std(faoyy5),N.std(faoyy6),N.std(faoyy7),N.std(faoyy8),N.std(faoyy9),N.std(faoyy10),N.std(faoyy11),N.std(faoyy12),N.std(faoyy13),N.std(faoyy14)]

faop=N.asarray(faop)/(10**6)
faop_std=N.asarray(faop_std)/(10**6)
allproisam=N.asarray(allproisam)/(10**6)
allproisam_std=N.asarray(allproisam_std)/(10**6)

#faoy=[2.7028,3.0312,2.8508,3.4374,3.4836,1.7038,2.1155,4.4007,3.3301,3.0681,2.6128,4.2432,3.0606,3.0640]
#faoa=[1560044,2376600,44712000,832000,10801214,26000,1903159,11793000,6302175,4038085,9891200,7666300,719370,698700]


 
fig = plt.figure(figsize=(20,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.4
opacity = 1.0
rects1 = plt.bar(index+0.1, allyisam, bar_width,yerr=allyisam_std,ecolor='black', capsize=5,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
#rects2 = plt.bar(index +0.18+ bar_width*1,allym3 , bar_width,
#                 alpha=opacity,
#                 color='g',
#                 label='M3')
rects0 = plt.bar(index + 0.1+bar_width*1, faoy, bar_width,yerr=faoy_std,ecolor='black', capsize=5,
                 alpha=opacity,
                 color='r',
                 label='Other studies')
#rects3 = plt.bar(index + 0.18+bar_width*3, iizuy1, bar_width,yerr=iizuy1_std,ecolor='black', capsize=5,
#                 alpha=opacity,
#                 color='y',
#                 label='Iizumi')

plt.xlabel('Country',fontsize=35)
plt.ylabel('Yield (t / ha cropland)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.11, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
leg=plt.legend(loc=1,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric2000yield1_countcombineflood3.png')
plt.show()



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
#rects2 = plt.bar(index +0.18+ bar_width*1,allprom3 , bar_width,
#                 alpha=opacity,
#                 color='g',
#                 label='M3')
rects0 = plt.bar(index + 0.1+bar_width*1, faop, bar_width,yerr=faop_std,ecolor='black', capsize=5,
                 alpha=opacity,
                 color='r',
                 label='Other studies')
#rects3 = plt.bar(index + 0.18+bar_width*3, iizup1, bar_width,yerr=iizup1_std,ecolor='black', capsize=5,
#                 alpha=opacity,
#                 color='y',
#                 label='iizumi')

plt.xlabel('Country',fontsize=35)
plt.ylabel('Production (Mtons)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.11, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
#leg=plt.legend(loc=2,ncol=3, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric2000pro_countluh2combineflood3.png')
plt.show()


fig = plt.figure(figsize=(20,10))
n_groups = 14

#ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.4
opacity = 1.0
rects1 = plt.bar(index+0.1, isamgrid, bar_width,yerr=isamgrid_std,ecolor='black', capsize=5,
                 alpha=opacity,
                 color='blue',
                 label='ISAM')
#rects2 = plt.bar(index +0.18+ bar_width*1,megrid , bar_width,
#                 alpha=opacity,
#                 color='g',
#                 label='M3')
rects0 = plt.bar(index + 0.1+bar_width*1, faoyy, bar_width,yerr=faoyy_std,ecolor='black', capsize=5,
                 alpha=opacity,
                 color='r',
                 label='Other studies')
#rects0 = plt.bar(index + 0.18+bar_width*3, iizuyg1, bar_width,yerr=iizuyg1_std,ecolor='black', capsize=5,
#                 alpha=opacity,
#                 color='y',
#                 label='Iizumi')

plt.xlabel('Country',fontsize=35)
plt.ylabel('(t grains / ha country)',fontsize=35)
#plt.ylabel('Yield (t / ha country)',fontsize=35)
plt.xticks(index + bar_width+0.11, ("Nepal","Pakistan","India","Sri Lanka","Bangladesh","Bhutan","Cambodia","Indonesia","Myanmar","Phillipines","Thailand","Vietnam","Laos","Malaysia"),rotation='vertical')
#lt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
#leg=plt.legend(loc=2,ncol=2, fontsize=35)
#leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=35)


plt.tight_layout()
plt.savefig('ric2000yield_countluh2combineflood3.png')
plt.show()

