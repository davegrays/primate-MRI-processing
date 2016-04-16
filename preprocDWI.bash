#!/bin/bash -e

if [ $# -lt 2 ]; then
	echo -e "\nUsage:	`basename $0` <groupFolder> <subjectName>"
	echo -e "will preprocess DWI of subject\n"
	echo -e "i.e	`basename $0` LES6 a_412\n"
	exit 1
fi

group=$1
sub=$2
maskfile="/media/data1/AMA09/Atlases/${group}/${sub}_mask"
rawdir="/media/data1/AMA09/Raw/${sub}/006__1009-1_2iso-NoG-mz-ep2d-minTE_DiCo/" #last '/' required

pushd ${group}/${sub}

case artifact_and_split in

artifact_and_split)
echo -e "\n$sub - finding the zero-artifact mask"
echo "spliting DWI"
fslsplit ${sub}_DWI
echo "reconstituting B0 and weighted series"
fslmerge -t ${sub}_B0s $(printf "vol%04d.nii.gz " $(seq -s " " 0 11))
fslmerge -t ${sub}_weighteds $(printf "vol%04d.nii.gz " $(seq -s " " 12 139))
rm vol*nii.gz
echo "finding max 10 image from weighted series"
fslmaths ${sub}_weighteds -Tmax -thr 10 -bin ${sub}_DWI_zeromask_dil

;&
upsamp)
echo -e "\nupsampling ${sub}_DWI.nii.gz to .6mm isotropic with spline interp. eliminate values < 0"
flirt -in ${sub}_DWI -ref ${sub}_DWI -applyisoxfm .6 -o DWI_upsamp -interp spline
fslmaths DWI_upsamp -thr 0 DWI_upsamp
echo "upsampling ${sub}_DWI_zeromask.nii.gz to .6mm isotropic linearly"
flirt -in ${sub}_DWI_zeromask -ref ${sub}_DWI_zeromask -applyisoxfm .6 -o DWI_zeromask
flirt -in ${sub}_DWI_zeromask_dil -ref ${sub}_DWI_zeromask_dil -applyisoxfm .6 -o DWI_zeromask_dil
fslmaths DWI_zeromask_dil -thr .5 -bin DWI_zeromask_dil
fslmaths DWI_upsamp -Tmin -bin -mas DWI_zeromask_dil DWI_zeromask

;&
coregister)
echo -e "\n$sub - putting T1 and brain mask in DWI native space" 
flirt -in ${sub}_T1 -ref DWI_upsamp -usesqform -applyxfm -omat T1_to_DWI.mat -interp spline -o T1_on_DWI
flirt -in ${maskfile} -ref DWI_upsamp -applyxfm -init T1_to_DWI.mat -o DWI_brainmask
fslmaths DWI_brainmask -thr .5 -bin DWI_brainmask

;&
bedpostX)
echo -e "\n$sub - getting mcverter bvecs. these are correct, while dcm2nii is Y-Z swapped"
rm -rf mcverter_out
mcverter -o mcverter_out -f fsl -d -n ${rawdir}
rm mcverter_out/*nii mcverter_out/*txt
echo -e "\n$sub - setting up bedpostX folder"
rm -rf bedpostX bedpostX.bedpostX
mkdir -p bedpostX
echo "converting bvals and bvecs (including the first 12 B0s)"
sed 's/^0/0 0 0 0 0 0 0 0 0 0 0 0/' mcverter_out/*bvecs > bedpostX/bvecs
sed 's/^0/0 0 0 0 0 0 0 0 0 0 0 0/' mcverter_out/*bvals > bedpostX/bvals
echo "make nodif_brain_mask - brainmask minus artifact voxels"
fslmaths DWI_brainmask -mas DWI_zeromask bedpostX/nodif_brain_mask
echo "link to DWI_upsamp"
ln -s `readlink -f DWI_upsamp.nii.gz` bedpostX/data.nii.gz
echo "running bedpostx"
bedpostx bedpostX --nf=2 --fudge=1 --bi=1000

;&
dtifit) ### THIS HASNT BEEN RUN YET (12/21/2014)
echo -e "\n$sub - getting diffusion tensor data"
mkdir -p dtifit
dtifit -k bedpostX/data.nii.gz -o dtifit/tensor -m bedpostX/nodif_brain_mask.nii.gz -r bedpostX/bvecs -b bedpostX/bvals

esac
popd
exit
