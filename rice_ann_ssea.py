from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import scipy.stats
from matplotlib.dates import DateFormatter
import datetime


country=NetCDFFile('/global/project/projectdirs/m1602/datasets4.full/arbit_init_state_05x05.nc','r')
#print iizumi
coun1 = country.variables['REGION_MASK'][:,:]
lon = country.variables['lon'][:]
coun,lona11=shiftgrid(180.5,coun1,lon,start=False)
print lona11


def annualyield(year,couna,counb):
    bb=year-1901
    aa=year-1901
    region1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/isamhiscru_rice_aogsluh2flood3.nc','r')
    maitrop = region1.variables['ricearea'][bb,:,:]
    maitrop=ma.masked_where(maitrop<=0,maitrop)
    maitrop=ma.filled(maitrop, fill_value=0.)

    maizeto = maitrop
    cc=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/rice_cheyenne/isamhiscru_rice_aogsluh2flood3.nc','r')
    isamcrumai3=cc.variables['yield'][aa,:,:]
    isamcrumai3=ma.masked_where(isamcrumai3<=0,isamcrumai3)
    isamcrumai3=ma.filled(isamcrumai3, fill_value=0.)
 
    yieldisamc3=0.

    a=0
    harea=0
    yieldic3=0.

    #USA 11501~11550
    for xx in range(0,360):
        for yy in range(0,720):
            if coun[xx,yy] >=couna and  coun[xx,yy] <=counb:
                harea=maizeto[xx,yy]+ harea
                yieldisamc3=(isamcrumai3[xx,yy]*maizeto[xx,yy])+yieldisamc3

                a=a+1
    yieldic3=N.average(isamcrumai3,weights=maizeto)

    return harea,yieldic3
    #return harea
def runmean(x,input):
        import pandas as pd
        #mean_zumiy1 = pd.rolling_mean(zumiy, window=5).shift(-2)
        meanout=pd.rolling_mean(input, window=5, center=True)
        mean_zumiy1=pd.rolling_mean(input, window=3, center=True)
        #print mean_zumiy1
        #print meanout
        meanout[1]=mean_zumiy1[1]
        meanout[x-2]=mean_zumiy1[x-2]
        meanout1=input-meanout
        return meanout1

#illzmui only 1983~2005
a1=1961
a2=2016
x=a2-a1
toarea= N.zeros(x)
zumiy= N.zeros(x)
clmy= N.zeros(x)
clmyn=N.zeros(x)
clmpn=N.zeros(x)
isamy= N.zeros(x)
zumip= N.zeros(x)
clmp= N.zeros(x)
isamp= N.zeros(x)
isamy1= N.zeros(x)
isamp1= N.zeros(x)
isamy2= N.zeros(x)
isamp2= N.zeros(x)
isamy3= N.zeros(x)
isamp3= N.zeros(x)
faoy=N.zeros(x)
faop=N.zeros(x)
#global!1981-2007
#faoy1=[34933.0,36090.0,29452.0,35256.0,37202.0,36279.,34863.,31001.,36186.,36907.,36996.,39064.,36305.,41120.,38100.,42061.,41511.,44354.,44255.,43236.,44775.,43884.,44610.,49452.,48196.,47719.,49980.]
#faop1=[446772517,448932280,347082034,450449992,485527301,478176622,453115794,403050234,476874503,483620724,494393020,533774898,477207493,568650520,517286851,586134845,584401847,615072804,607426254,592030667,6155143531,603544019,645048171,729511789,714185792,707932497,793055503]

#1961-2015
#faoy1=[18693,18958,20567,21025,20353,20781,21749,22328,22540,23808,23619,23245,24524,24251,25186,24509,25718,26843,26591,27482,28272,29804,31367,32261,32570,32441,32651,33296,34541,35286,35362,35847,36162,36576,36580,37840,38165,38160,38970,38874,39502,38630,39538,40301,40848,41187,42273,42920,43453,43364,44640,45396,45097,45572,46036]
faoy1=[1.595448629,1.524795006,1.641391236,1.660123143,1.503442743,1.482760268,1.624998526,1.710303606,1.753404918,1.834178929,1.834934429,1.762638566,1.883829891,1.804949704,1.950771346,1.878620677,2.004723,2.073161917,1.92211414,2.181764345,2.217288568,2.227052277,2.418344642,2.420286446,2.523002725,2.482690932,2.48643258,2.661728653,2.760806716,2.772124325,2.806339672,2.829596894,2.946157054,2.975380228,2.937411503,3.03031447,3.034448683,3.06730319,3.171623559,3.208767222,3.337266602,3.20306261,3.431566164,3.399451518,3.531234278,3.55763943,3.673805641,3.70941846,3.740051595,3.806387299,3.920808979,4.014002195,4.002307089,4.011326448,4.044533691]
faoy=N.asarray(faoy1)

#us
#faoy1=[68380,71082,50893,66981,74074,74930,75225,53109,72978,74380,68172,82526,63211,86997,71230,79777,79522,84382,83979,85910,86733,81179,89246,100636,92852]
#faop1=[206222000,209180000,106030000,194880000,225453008,208943008,181142000,125194000,191319008,201532000,189866496,240719008,160984992,255292992,187968992,234527008,233867008,247882000,239548992,251852210,241375035,227765357,256227304,299873563,282260662]


#name=["Global","USA","China","Brazil","Argentina","Mexico","India","Ukraine","Indonesia","France","Southafrica"]
#range1=[10100,11501,41501,20301,20101,11101,42901,47800,50300,42200,33900]
#range2=[50700,11550,41529,20327,20124,11132,42929,47800,50300,42200,33900]

#name=["Italy","Canada","Vietnam","Hungary","Romania","Philippines","Thailand","Chile","Spain","Nigeria","Germany"]
#range1=[43400,10201,48200,42700,46000,50700,47500,20400,46800,33400,42500]
#range2=[43400,10212,48200,42700,46000,50700,47500,20400,46800,33400,42500]

name=["SSEA"]
range1=[7.9]
range2=[8.1]
#name=["US"]
#range1=[11501]
#range2=[11550]

for i, name1 in enumerate(name):
        a=0
        clmy1=N.zeros(x)
	toarea= N.zeros(x)
	zumiy= N.zeros(x)
	clmy= N.zeros(x)
	isamy= N.zeros(x)
	zumip= N.zeros(x)
	clmp= N.zeros(x)
	isamp= N.zeros(x)
	isamy1= N.zeros(x)
	isamp1= N.zeros(x)
        isamy2= N.zeros(x)
        isamp2= N.zeros(x)
        isamy3= N.zeros(x)
        isamp3= N.zeros(x)

        clmyn=N.zeros(x)
        clmpn=N.zeros(x)
	for num in range(a1,a2):
	    reu=annualyield(num,range1[i],range2[i])
	    toarea[a]=reu[0]
            isamy3[a]=reu[1]

	    a=a+1

        aisamy3=isamy3-N.average(isamy3)

        afaoy=faoy-N.average(faoy)

        mean_isamy3=runmean(x,isamy3)

        mean_faoy=runmean(x,faoy)



	fig = plt.figure(figsize=(26,20))


	ax = fig.add_subplot(321)
	xx=range(a1,a2)
	x1=range(1902,2016)
	x2=range(1982,2007)
	xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
        xdates2 = [datetime.datetime.strptime(str(int(date)),'%Y') for date in x2]

	ax.plot_date(xdates,faoy,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,isamy3,"ko-",label="ISAM")


	ax.xaxis.set_major_formatter(DateFormatter('%Y'))
        ccp3=scipy.stats.pearsonr(isamy3,faoy)

        leg=plt.legend(['FAO','ISAM {:04.2f}'.format(ccp3[0])],fontsize=18,loc=4)

#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)

	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Rice yield (t/ha)",fontsize=18)
	plt.title("Yield",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)





	ax = fig.add_subplot(325)
        ax.plot_date(xdates,mean_faoy,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,mean_isamy3,"ko-",label="ISAM")
        ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	mean1_faoy=N.zeros([x-2])
        mean1_isamy3=N.zeros([x-2])

	for i in range(1,x-2):
		mean1_faoy[i]=mean_faoy[i]
                mean1_isamy3[i]=mean_isamy3[i]

        ccp4=scipy.stats.pearsonr(mean1_isamy3,mean1_faoy)

        leg=plt.legend(['FAO', 'ISAM {:04.2f}'.format(ccp4[0])],fontsize=18,loc=4)

	leg.get_frame().set_alpha(0.5)


#	plt.title("Anormaly by subtracting a moving mean average",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Rice yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)

        ax = fig.add_subplot(326)
        ax.plot_date(xdates2,mean_faoy[21:46],"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates2,mean_isamy3[21:46],"ko-",label="ISAM")
        ax.xaxis.set_major_formatter(DateFormatter('%Y'))



        ccp4=scipy.stats.pearsonr(mean_faoy[21:46],mean_isamy3[21:46])

        leg=plt.legend(['FAO', 'ISAM {:04.2f}'.format(ccp4[0])],fontsize=18,loc=4)

        leg.get_frame().set_alpha(0.5)


#       plt.title("Anormaly by subtracting a moving mean average",fontsize=18)
        plt.xlabel("Year",fontsize=18)
        plt.ylabel("Rice yield (t/ha)",fontsize=18)
        plt.tick_params(axis='both',labelsize=18)



	plt.savefig('rice_fao_{0}_1961_2016.png'.format(name1),dpi=600,bbox_inches='tight')
plt.show()


