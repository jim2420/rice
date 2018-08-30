#!/bin/sh
module load cdo

echo "writing out yield"



name=(equilibrium equilibrium_cli equilibrium_co2 equilibrium_fert equilibrium_irr equilibrium_ndep ric_cli ric_co2 ric_fert ric_irr ric_irr_fert ric_ndep)
for i in {0..11}
do
 for f in {1901..2015}
  do
   echo $f ${name[i]}  
#  cd /scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/${name[i]}/output
  cd /scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/experiments/${name[i]}/output
   echo *_crop_$f.nc crop_$f.nc
#   cdo selname,totalyield *_crop_$f.nc crop_$f.nc
   cdo selname,ch4_flux,ch4_oxid *.bgc-yearly-2d_$f.nc crop_$f.nc
   cdo setyear,$f crop_$f.nc sam_$f.nc
  done
cdo mergetime sam_*.nc ${name[i]}_ch4.nc
rm sam_*.nc
rm crop_*.nc
done

 
