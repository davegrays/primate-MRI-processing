#!/bin/bash -e

if [ $# -lt 5 ]; then
	echo -e "\nUsage:	`basename $0` <ageFolder> <file_list> <initial_target_HEAD> <initial_target_BRAIN> <num_cores>"
	echo -e "will make template using parallelization on num_cores\n"
	echo -e "i.e	`basename $0` 260w dicom_to_nifti_3T_lists/3T_5yr_list_niftis.txt MASTER_template/template_HEAD.nii.gz MASTER_template/template_BRAIN.nii.gz 6\n"
	exit 1
fi

age=$1
sublist=`cat $2 | tr "\n" " "`
targ=`readlink -f $3`
targ_brain=`readlink -f $4`
num_cores=$5

pushd $age

#use dof6 to first get the brain masks
#use dof6_2 if you already have the brain masks

case dof6_2 in

dof6)
#put images in template coordinate space using whole head
#using flirt because ANTS takes forever on this stage for some reason
parallel -j ${num_cores} -r "flirt -dof 6 -in * -ref ${targ} -o *_rot2atl -omat *_rot2atl.mat -interp spline" ${sublist}

;&
dof12)
#create linear affine template (no bias correction)
buildtemplateparallel.sh -d 3 -n 0 -i 1 -m 1x0x0 -c 2 -j ${num_cores} -o ${age}_affine *rot2atl.nii.gz
mkdir -p affine_avg
mv ${age}* *cfg affine_avg
rm affine_avg/*Warp.nii.gz affine_avg/*repaired.nii.gz
rm -r GR_iteration_0

;&
full)
#create deformable template, only 1 iteration here
buildtemplateparallel.sh -d 3 -z affine_avg/${age}_affinetemplate.nii.gz -i 1 -m 50x90x30 -c 2 -j ${num_cores} -o ${age}HEAD *rot2atl.nii.gz
mv *Affine.txt GR_iteration_0
rm ${age}* *cfg

;&
get_brain_masks)
#target mask
fslmaths ${targ_brain} -bin templateMask
#get warp from target to template
ANTS 3 -m  CC[GR_iteration_0/${age}HEADtemplate.nii.gz,${targ},1,5] -t SyN[0.25] -r Gauss[3,0] -o targ_to_${age}HEADtemplate -i 50x90x30 --use-Histogram-Matching  --number-of-affine-iterations 10000x10000x10000x10000x10000 --MI-option 32x16000
#apply warp from target to template to each subject
parallel -j ${num_cores} -r "antsApplyTransforms -d 3 -i templateMask.nii.gz -t [GR_iteration_0/${age}HEAD*_rot2atlAffine.txt,1] GR_iteration_0/${age}HEAD*_rot2atlInverseWarp.nii.gz targ_to_${age}HEADtemplateWarp.nii.gz targ_to_${age}HEADtemplateAffine.txt -r *_rot2atl.nii.gz -o *_rot2atl_mask.nii.gz" ${sublist}
#back it into native space; threshold it; apply it
parallel -j ${num_cores} -r "convert_xfm -omat *_rot2native.mat -inverse *_rot2atl.mat; flirt -in *_rot2atl_mask -ref * -o *_mask -applyxfm -init *_rot2native.mat; fslmaths *_mask -thr .5 -bin *_mask; fslmaths * -mas *_mask *_brain" ${sublist}
rm *_rot2native.mat *_rot2atl_mask.nii.gz templateMask.nii.gz

;&
dof6_2)
#using brain and applying to nonbrain
parallel -j ${num_cores} -r "flirt -dof 6 -in *_brain -ref ${targ_brain} -omat *_rot2atl.mat" ${sublist}
parallel -j ${num_cores} -r "flirt -in * -ref ${targ} -o *_rot2atl -applyxfm -init *_rot2atl.mat -interp spline" ${sublist}

;&
dof12_2)
rm -r affine_avg
buildtemplateparallel.sh -d 3 -n 0 -i 1 -m 1x0x0 -c 2 -j ${num_cores} -o ${age}_affine *rot2atl.nii.gz
mkdir -p affine_avg
mv ${age}* *cfg affine_avg
rm affine_avg/*Warp.nii.gz affine_avg/*repaired.nii.gz
rm -r GR_iteration_0

;&
full_2)
buildtemplateparallel.sh -d 3 -z affine_avg/${age}_affinetemplate.nii.gz -m 50x90x30 -c 2 -j ${num_cores} -o ${age}HEAD *rot2atl.nii.gz

;&
get_brain_masks_2)
#target mask
fslmaths ${targ_brain} -bin templateMask
#get warp from target to template
ANTS 3 -m  CC[${age}HEADtemplate.nii.gz,${targ},1,5] -t SyN[0.25] -r Gauss[3,0] -o targ_to_${age}HEADtemplate -i 50x90x30 --use-Histogram-Matching  --number-of-affine-iterations 10000x10000x10000x10000x10000 --MI-option 32x16000
#apply warp from target to template to each subject
parallel -j ${num_cores} -r "antsApplyTransforms -d 3 -i templateMask.nii.gz -t [${age}HEAD*_rot2atlAffine.txt,1] ${age}HEAD*_rot2atlInverseWarp.nii.gz targ_to_${age}HEADtemplateWarp.nii.gz targ_to_${age}HEADtemplateAffine.txt -r *_rot2atl.nii.gz -o *_rot2atl_mask.nii.gz" ${sublist}
#back it into native space; threshold it; apply it
parallel -j ${num_cores} -r "convert_xfm -omat *_rot2native.mat -inverse *_rot2atl.mat; flirt -in *_rot2atl_mask -ref * -o *_mask -applyxfm -init *_rot2native.mat; fslmaths *_mask -thr .5 -bin *_mask; fslmaths * -mas *_mask *_brain" ${sublist}
rm *_rot2native.mat *_rot2atl_mask.nii.gz templateMask.nii.gz

esac

popd

exit
