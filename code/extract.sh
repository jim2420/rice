#!/bin/sh
module load cdo

echo "writing out yield"



name=(ric_irr_fert_long4 ric_fertall ric_fert ric_irr_fert)
for i in {2..3}
do
 for f in {1901..2015}
  do
   echo $f ${name[i]}  
#  cd /scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/${name[i]}/output
  cd /scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/india/${name[i]}/output
   echo *_crop_$f.nc crop_$f.nc
   cdo selname,totalyield *_crop_$f.nc crop_$f.nc
   cdo setyear,$f crop_$f.nc sam_$f.nc
  done
cdo mergetime sam_*.nc ${name[i]}.nc
rm sam_*.nc
rm crop_*.nc
done

 
