#!/bin/sh
module load cdo

echo "writing out yield"



name=(ric_fert ric_irr_fert)
for i in {0..1}
do
 for f in {1901..2015}
  do
   echo $f ${name[i]}  
#  cd /scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/${name[i]}/output
#  cd /scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/nonflood/${name[i]}/output
  cd /scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/flood3/${name[i]}/output

   echo *_crop_$f.nc crop_$f.nc
   cdo selname,ch4_flux,ch4_oxid *.bgc-yearly-2d_$f.nc crop_$f.nc
   cdo setyear,$f crop_$f.nc sam_$f.nc
  done
cdo mergetime sam_*.nc ${name[i]}_ch4.nc
rm sam_*.nc
rm crop_*.nc
done

 
