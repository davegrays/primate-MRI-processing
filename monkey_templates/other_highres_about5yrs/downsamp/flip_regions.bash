#!/bin/bash

## FOR THIS SCRIPT TO WORK YOU MUST START WITH HANDDRAWN REGIONS ON THE LEFT SIDE
## NAME THEM <ROI>_L.nii.gz
## THIS WILL OVERWRITE ALL THE <ROI>_R.nii.gz files 

inROIdir="ROIs_subcort"
refbrain="CON8HEADtemplate_BRAIN.nii.gz"

for ROIpath in ${inROIdir}/*_L.nii.gz;do
	#declare string name variables
	ROI_L_niigz=`echo ${ROIpath} | sed 's/.*\///'` #"ROI_L.nii.gz"
	ROI_R_niigz=`echo ${ROI_L_niigz} | sed 's/_L\.nii\.gz/_R\.nii\.gz/'` #"ROI_R.nii.gz"
	ROI_base=`echo ${ROI_L_niigz} | sed 's/_L\.nii\.gz//'` #"ROI"
	#flip the roi
	fslswapdim ${inROIdir}/${ROI_L_niigz} -x y z ${inROIdir}/${ROI_base}_tmpflip
	#warp the flipped roi to the orig image. using linear interp instead of nn
	antsApplyTransforms -d 3 -i ${inROIdir}/${ROI_base}_tmpflip.nii.gz -t ../asymmetry/flipped_to_origWarp.nii.gz ../asymmetry/flipped_to_origAffine.txt -r ${refbrain} -o ${inROIdir}/${ROI_R_niigz}
	#threshold at .5
	fslmaths ${inROIdir}/${ROI_R_niigz} -thr .5 -bin ${inROIdir}/${ROI_R_niigz}
	#clean up
	rm ${inROIdir}/${ROI_base}_tmpflip.nii.gz
done
exit
