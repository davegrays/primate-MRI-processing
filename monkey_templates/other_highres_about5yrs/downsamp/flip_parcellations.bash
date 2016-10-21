#!/bin/bash

## FOR THIS SCRIPT TO WORK YOU MUST START WITH MANUALLY EDITED SEGMENTATION/PARCELLATION FILES ON LEFT SIDE
## NAME THEM <ParcName>_L.nii.gz
## you should have files mask_hardsegL.nii.gz and mask_hardsegR.nii.gz for separating left and right segs/parcs
## N.B. in parcellation files, regional intensity values should be ordered 1, 2, 3, ..., maxVal (with no skips or repeats)
## N.B. parcellation files must have less than 10000 regions

inROIdir="flip_parcs"
refbrain="CON8HEADtemplate_BRAIN.nii.gz"
cp HardSeg.nii.gz ${inROIdir}/HardSeg_L.nii.gz
cp MY_Subcorts.nii.gz ${inROIdir}/MY_Subcorts_L.nii.gz

case phase1 in

phase1)
for ROIpath in ${inROIdir}/*_L.nii.gz;do
	#get max value of file (rounded to nearest integer)
	nt=`fslstats ${ROIpath} -R | awk '{print $2}' | awk '{printf("%d\n",$1 + 0.5)}'`

	#declare string name variables
	ROI_L_niigz=`echo ${ROIpath} | sed 's/.*\///'` #"ROI_L.nii.gz"
	ROI_R_niigz=`echo ${ROI_L_niigz} | sed 's/_L\.nii\.gz/_R\.nii\.gz/'` #"ROI_R.nii.gz"
	ROI_base=`echo ${ROI_L_niigz} | sed 's/_L\.nii\.gz//'` #"ROI"

	#flip the roi
	fslswapdim ${inROIdir}/${ROI_L_niigz} -x y z ${inROIdir}/${ROI_base}_tmpflip

	#binarize, project, and take 50%threshold (final mask)
	fslmaths ${inROIdir}/${ROI_base}_tmpflip -bin ${inROIdir}/${ROI_base}_flipmask
	antsApplyTransforms -d 3 -i ${inROIdir}/${ROI_base}_flipmask.nii.gz -t ../asymmetry/flipped_to_origWarp.nii.gz ../asymmetry/flipped_to_origAffine.txt -r ${refbrain} -o ${inROIdir}/${ROI_base}_mask_R.nii.gz
	fslmaths ${inROIdir}/${ROI_base}_mask_R.nii.gz -thr .5 -bin ${inROIdir}/${ROI_base}_mask_R.nii.gz

	#extract and warp each tissue class separately
	for ((c=1;c<=$nt;c++));do
		n=`printf %05d $c`
		fslmaths ${inROIdir}/${ROI_base}_tmpflip -thr $c -uthr $c -bin ${inROIdir}/${ROI_base}_tmpflip_$n
		#warp the flipped roi to the orig image. using linear interp instead of nn
		antsApplyTransforms -d 3 -i ${inROIdir}/${ROI_base}_tmpflip_$n.nii.gz -t ../asymmetry/flipped_to_origWarp.nii.gz ../asymmetry/flipped_to_origAffine.txt -r ${refbrain} -o ${inROIdir}/${ROI_base}_${n}_R.nii.gz
	done
	#recombine split files into 4D file
	fslmerge -t ${inROIdir}/${ROI_base}_R_4D ${inROIdir}/${ROI_base}_0*nii.gz
	#convert 4D file into ranks (and add 1 because fsl uses 0 starting index)
	fslmaths ${inROIdir}/${ROI_base}_R_4D -rank -Tmaxn -add 1 -mas ${inROIdir}/${ROI_base}_mask_R.nii.gz ${inROIdir}/${ROI_R_niigz}
	#clean up
	rm ${inROIdir}/${ROI_base}_tmpflip*nii.gz ${inROIdir}/${ROI_base}_0*nii.gz
	rm ${inROIdir}/${ROI_base}_flipmask.nii.gz ${inROIdir}/${ROI_base}_mask_R.nii.gz ${inROIdir}/${ROI_base}_R_4D.nii.gz

	############
	# RINSE AND REPEAT
	# project right side back onto left (allows for interpolation-related blurring)
	############

	#flip the roi
#	fslswapdim ${inROIdir}/${ROI_R_niigz} -x y z ${inROIdir}/${ROI_base}_tmpflip

	#binarize, project, and take 50%threshold (final mask)
#	fslmaths ${inROIdir}/${ROI_base}_tmpflip -bin ${inROIdir}/${ROI_base}_flipmask
#	antsApplyTransforms -d 3 -i ${inROIdir}/${ROI_base}_flipmask.nii.gz -t ../asymmetry/flipped_to_origWarp.nii.gz ../asymmetry/flipped_to_origAffine.txt -r ${refbrain} -o ${inROIdir}/${ROI_base}_mask_L.nii.gz
#	fslmaths ${inROIdir}/${ROI_base}_mask_L.nii.gz -thr .5 -bin ${inROIdir}/${ROI_base}_mask_L.nii.gz

	#extract and warp each tissue class separately
#	for ((c=1;c<=$nt;c++));do
#		n=`printf %05d $c`
#		fslmaths ${inROIdir}/${ROI_base}_tmpflip -thr $c -uthr $c -bin ${inROIdir}/${ROI_base}_tmpflip_$n
		#warp the flipped roi to the orig image. using linear interp instead of nn
#		antsApplyTransforms -d 3 -i ${inROIdir}/${ROI_base}_tmpflip_$n.nii.gz -t ../asymmetry/flipped_to_origWarp.nii.gz ../asymmetry/flipped_to_origAffine.txt -r ${refbrain} -o ${inROIdir}/${ROI_base}_${n}_L.nii.gz
#	done
	#recombine split files into 4D file
#	fslmerge -t ${inROIdir}/${ROI_base}_L_4D ${inROIdir}/${ROI_base}_0*nii.gz
	#convert 4D file into ranks (and add 1 because fsl uses 0 starting index)
#	fslmaths ${inROIdir}/${ROI_base}_L_4D -rank -Tmaxn -add 1 -mas ${inROIdir}/${ROI_base}_mask_L.nii.gz ${inROIdir}/${ROI_L_niigz}
	#clean up
#	rm ${inROIdir}/${ROI_base}_tmpflip*nii.gz ${inROIdir}/${ROI_base}_0*nii.gz
#	rm ${inROIdir}/${ROI_base}_flipmask.nii.gz ${inROIdir}/${ROI_base}_mask_L.nii.gz ${inROIdir}/${ROI_base}_L_4D.nii.gz

done

;&
phase2)
echo "masking left and right separately"
echo "fslmaths ${inROIdir}/HardSeg_L -mas ${inROIdir}/mask_hardsegL ${inROIdir}/HardSeg_L_masked"
fslmaths ${inROIdir}/HardSeg_L -mas ${inROIdir}/mask_hardsegL ${inROIdir}/HardSeg_L_masked
echo "fslmaths ${inROIdir}/HardSeg_R -mas ${inROIdir}/mask_hardsegR ${inROIdir}/HardSeg_R_masked"
fslmaths ${inROIdir}/HardSeg_R -mas ${inROIdir}/mask_hardsegR ${inROIdir}/HardSeg_R_masked
echo "fslmaths ${inROIdir}/MY_Subcorts_L -mas ${inROIdir}/mask_MY_SubcortsL ${inROIdir}/subcorts_L_masked"
fslmaths ${inROIdir}/MY_Subcorts_L -mas ${inROIdir}/mask_hardsegL ${inROIdir}/subcorts_L_masked
echo "fslmaths ${inROIdir}/MY_Subcorts_R -mas ${inROIdir}/mask_MY_SubcortsR ${inROIdir}/subcorts_R_masked"
fslmaths ${inROIdir}/MY_Subcorts_R -mas ${inROIdir}/mask_hardsegR ${inROIdir}/subcorts_R_masked

echo "recombining left and right hems into one"
echo "fslmaths ${inROIdir}/HardSeg_R_masked -add ${inROIdir}/HardSeg_L_masked HardSeg_BOTH"
fslmaths ${inROIdir}/HardSeg_R_masked -add ${inROIdir}/HardSeg_L_masked ${inROIdir}/HardSeg_BOTH
echo "fslmaths ${inROIdir}/subcorts_R_masked -add ${inROIdir}/subcorts_L_masked MY_Subcorts_BOTH"
fslmaths ${inROIdir}/subcorts_R_masked -add ${inROIdir}/subcorts_L_masked ${inROIdir}/MY_Subcorts_BOTH

;&
phase3)
echo "copying to FINAL_ROIS folder"
cp ${inROIdir}/MY_Subcorts_BOTH.nii.gz FINAL_ROIS/MY_Subcorts_BOTH.nii.gz
cp ${inROIdir}/HardSeg_BOTH.nii.gz FINAL_ROIS/HardSeg_BOTH.nii.gz

echo "executing extract_regions.bash and extract_tissues.bash"
pushd FINAL_ROIS
./extract_regions.bash
./extract_tissues.bash
popd

esac

exit
