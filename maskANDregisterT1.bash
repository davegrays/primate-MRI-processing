#!/bin/bash -e

if [ $# -ne 4 ]; then
	echo -e "\nUsage:	`basename $0` <subject_directory> <subject_T1_prefix> <atlas_HEAD_prefix> <atlas_BRAIN_prefix>"
	echo -e "i.e	`basename $0`  /path/to/monkeyAdir monkeyA /path/to/master_template_HEAD /path/to/master_template_BRAIN\n"
	echo -e "will process the file in /path/to/monkeyAdir/monkeyA.nii.gz"
	echo -e "atlas files should also have .nii.gz format\n"
	exit 1
fi

dir=$1
sub=$2
atl=`readlink -f $3`
atl_brain=`readlink -f $4`

#check for existence of c3d_affine_tool in path
path2c3d=`which c3d_affine_tool`
if [[ $path2c3d == "" ]]; then
	echo "FATAL ERROR: c3d_affine_tool not found."
	echo "Make sure convert3D is properly installed and added to your path before using this tool."
	exit 1
fi

pushd $dir

#use dof6_head to first get the brain masks
#use dof6_brain if you already have the brain masks

case dof6_head in

dof6_head)
#put images in template coordinate space using whole head
#using flirt because ANTS takes forever on this stage for some reason
flirt -v -dof 6 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -in ${sub} -ref ${atl} -o ${sub}_rot2atl -omat ${sub}_rot2atl.mat -interp sinc

;&
get_brain_masks)
#atlas mask
fslmaths ${atl_brain} -bin templateMask
#get warp from atlas head to subject head
ANTS 3 -m  CC[${sub}_rot2atl.nii.gz,${atl}.nii.gz,1,5] -t SyN[0.25] -r Gauss[3,0] -o atl2T1rot -i 60x50x20 --use-Histogram-Matching  --number-of-affine-iterations 10000x10000x10000x10000x10000 --MI-option 32x16000
#apply warp
antsApplyTransforms -d 3 -i templateMask.nii.gz -t atl2T1rotWarp.nii.gz atl2T1rotAffine.txt -r ${sub}_rot2atl.nii.gz -o ${sub}_rot2atl_mask.nii.gz
## NOTE: antsApplyTransforms is deprecated in ANTS 2.0 and might not be available; if so, use WarpImageMultiTransform instead
#back it into native space; threshold it; apply it
convert_xfm -omat ${sub}_rot2native.mat -inverse ${sub}_rot2atl.mat; flirt -in ${sub}_rot2atl_mask -ref ${sub} -o ${sub}_mask -applyxfm -init ${sub}_rot2native.mat
fslmaths ${sub}_mask -thr .5 -bin ${sub}_mask
fslmaths ${sub} -mas ${sub}_mask ${sub}_brain
rm ${sub}_rot2native.mat ${sub}_rot2atl_mask.nii.gz templateMask.nii.gz

;&
dof6_brain)
flirt -v -dof 6 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -in ${sub}_brain -ref ${atl_brain} -o ${sub}_brain_rot2atl -omat ${sub}_rot2atl.mat
flirt -in ${sub} -ref ${atl} -o ${sub}_rot2atl -applyxfm -init ${sub}_rot2atl.mat -interp sinc
#convert FSL mat to ANTS compatible itk mat
c3d_affine_tool -ref ${sub}_brain_rot2atl.nii.gz -src ${sub}_brain.nii.gz ${sub}_rot2atl.mat -fsl2ras -oitk ${sub}_rot2atl.txt
#N.B. if you want to use the inverse, specify so in the antsApplyTransforms command (rather than creating the inverse in FSL and converting to itk, which works improperly)

;&
final_warp)
#get warp from atlas brain to subject brain
ANTS 3 -m  CC[${sub}_brain_rot2atl.nii.gz,${atl_brain}.nii.gz,1,5] -t SyN[0.25] -r Gauss[3,0] -o atl2T1rot -i 50x60x30 --use-Histogram-Matching  --number-of-affine-iterations 10000x10000x10000x10000x10000 --MI-option 32x16000
#apply warp
antsApplyTransforms -d 3 -i ${atl_brain}.nii.gz -t atl2T1rotWarp.nii.gz atl2T1rotAffine.txt -r ${sub}_rot2atl.nii.gz -o atl2T1rot_deforemdImage.nii.gz
## NOTE: antsApplyTransforms is deprecated in ANTS versions > 2.0 and might not be available; if so, use WarpImageMultiTransform instead

esac

popd

exit
