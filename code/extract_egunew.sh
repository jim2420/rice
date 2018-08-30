#!/bin/sh
module load cdo

echo "writing out yield"


#name=(equilibrium_irr_co2 equilibrium_irr_fert equilibrium_nondep equilibrium_cli equilibrium_fert equilibrium equilibrium_co2 equilibrium_irr)
#name=(ric_fert  ric_fert_constcli  ric_irr  ric_irr_fert  ric_irr_fert_constc  ric_irr_fert_constcli  ric_irrfixed_fert)
name=(ric_fert ric_irr_fert)
for i in {0..1}
do
 for f in {1901..2015}
  do
   echo $f ${name[i]}  
  cd /scratch2/scratchdirs/tslin2/isam/rice_cheyenne/his_cru/aogs/flood1/${name[i]}/output
   echo *_crop_$f.nc crop_$f.nc
   cdo selname,totalyield *_crop_$f.nc crop_$f.nc
   cdo setyear,$f crop_$f.nc sam_$f.nc
  done
cdo mergetime sam_*.nc ${name[i]}.nc
rm sam_*.nc
rm crop_*.nc
done

 
