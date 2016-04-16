#!/bin/bash -e

if [ $# -lt 2 ]; then
	echo -e "\nUsage:	`basename $0` <groupFolder> <subjectName>"
	echo -e "will generate atlas-based segmentations and parcellations in DWI space of subject. then runs probtrackX for connectomes\n"
	echo -e "i.e	`basename $0` LES6 a_412\n"
	exit 1
fi

group=$1
sub=$2
rotdir="/media/data1/AMA09/Atlases/${group}"
warpdir="/media/data1/AMA09/Atlases/${group}/subjectwarps"
atldir="/media/data1/AMA09/Atlases/CON8/downsamp/FINAL_ROIS"
cwd=$PWD

pushd ${group}/${sub}
mkdir -p segmentation parcellation

case prep_probtrack in

segment)
c3d_affine_tool -ref ${rotdir}/${sub}_rot2atl.nii.gz -src ${rotdir}/${sub}.nii.gz ${rotdir}/${sub}_rot2atl.mat -fsl2ras -oitk ${sub}_T1rot2atl.txt
c3d_affine_tool -ref T1_on_DWI.nii.gz -src ${sub}_T1.nii.gz T1_to_DWI.mat -fsl2ras -oitk T1_to_DWI.txt

echo -e "\n$sub - propagating tissues to DWI native space"
c=0
for tiss in CSF GMcort WM GMdeep;do
	c=$(($c+1))
	n=`printf %05d $c`
	antsApplyTransforms -d 3 -i ${atldir}/tissue_${tiss}.nii.gz -t T1_to_DWI.txt [${sub}_T1rot2atl.txt,1] ${warpdir}/${sub}_rot2atlWarp.nii.gz ${warpdir}/${sub}_rot2atlAffine.txt -r T1_on_DWI.nii.gz -o segmentation/tissue_${n}.nii.gz
done
antsApplyTransforms -d 3 -i ${atldir}/HardSeg_mask.nii.gz -t T1_to_DWI.txt [${sub}_T1rot2atl.txt,1] ${warpdir}/${sub}_rot2atlWarp.nii.gz ${warpdir}/${sub}_rot2atlAffine.txt -r T1_on_DWI.nii.gz -o segmentation/hardseg_mask.nii.gz
fslmaths segmentation/hardseg_mask.nii.gz -thr .5 -bin segmentation/hardseg_mask.nii.gz

echo "combining tissue probabilities into hardseg"
#recombine split files into 4D file
fslmerge -t segmentation/tissue_4D segmentation/tissue_0*nii.gz
#convert 4D file into ranks (and add 1 because fsl uses 0 starting index)
fslmaths segmentation/tissue_4D -rank -Tmaxn -add 1 -mas segmentation/hardseg_mask.nii.gz segmentation/hardseg
#clean up
rm segmentation/tissue_0*nii.gz segmentation/tissue_4D.nii.gz

echo "re-split into individual binarized tissues"
c=0
for tiss in CSF GMcort WM GMdeep;do
	c=$(($c+1))
	fslmaths segmentation/hardseg -thr $c -uthr $c -bin segmentation/${tiss}
done

;&
make_cortParc)
#you go from F99 to CON8 to subject
inparc="/media/data1/Primate_Atlases/RM/RM_F99_ROItemplate.nii.gz"
refdir="/media/data1/AMA09/Atlases/CON8"

#f99 to CON8 to subject
antsApplyTransforms -d 3 -i ${inparc} -t T1_to_DWI.txt [${sub}_T1rot2atl.txt,1] ${warpdir}/${sub}_rot2atlWarp.nii.gz ${warpdir}/${sub}_rot2atlAffine.txt ${refdir}/F99brain_to_CON8Warp.nii.gz ${refdir}/F99brain_to_CON8Affine.txt -r T1_on_DWI.nii.gz -o parcellation/RMorig.nii.gz -n NearestNeighbor

#mask the parcellation on the cortGM tissue class and dilate 5 times, preserving the boundaries of each previous iteration
echo "masking and dilating cortGM parcellation"
fslmaths parcellation/RMorig -mas segmentation/GMcort parcellation/RMorig
echo "dilate 1"
fslmaths parcellation/RMorig -dilD -mas segmentation/GMcort parcellation/RMdil1
fslmaths parcellation/RMorig -bin -mul -10000 -add parcellation/RMdil1 -thr 0 -add parcellation/RMorig parcellation/RMdil1
for c in 2 3 4 5;do
	echo "dilate $c"
	p=$(($c-1))
	fslmaths parcellation/RMdil$p -dilD -mas segmentation/GMcort parcellation/RMdil$c
	fslmaths parcellation/RMdil$p -bin -mul -10000 -add parcellation/RMdil$c -thr 0 -add parcellation/RMdil$p parcellation/RMdil$c
done
mv parcellation/RMdil${c}.nii.gz parcellation/RMfinal.nii.gz
rm parcellation/RMdil*.nii.gz

;&
split_corticals)
#split the ROIs into individual files according to LUT file
LUT="/media/data1/Primate_Atlases/RM/MYedited_RegionNames.txt"
declare -a val=($(cat ${LUT} | sed "1d" | awk '{print $1}' | tr "\n" " "))
declare -a reg=($(cat ${LUT} | sed "1d" | awk '{print $2}' | tr "\n" " "))
end=${#val[@]}
end=$(($end - 1)) #because of 0 indexing
mkdir -p parcellation/ROIS_full
for ((n=0;n<=$end;n++));do
	echo "fslmaths parcellation/RMfinal -thr ${val[$n]} -uthr ${val[$n]} -bin parcellation/ROIS_full/${reg[$n]}"
	fslmaths parcellation/RMfinal -thr ${val[$n]} -uthr ${val[$n]} -bin parcellation/ROIS_full/${reg[$n]}
done

#do the same for the inner cortGM edge
fslmaths segmentation/WM -dilD segmentation/WM_dil
fslmaths parcellation/RMfinal -mas segmentation/WM_dil parcellation/RMfinal_edge
mkdir -p parcellation/ROIS_edge
for ((n=0;n<=$end;n++));do
	echo "fslmaths parcellation/RMfinal_edge -thr ${val[$n]} -uthr ${val[$n]} -bin parcellation/ROIS_edge/${reg[$n]}"
	fslmaths parcellation/RMfinal_edge -thr ${val[$n]} -uthr ${val[$n]} -bin parcellation/ROIS_edge/${reg[$n]}
done

;&
parcellate_subcorts)
mkdir -p parcellation/subcorts_full parcellation/subcorts_edge
echo -e "\n$sub - propagating subcortical regions to DWI native space"
c=0
for tiss in amyg_R caud_R GP_R hip_R NAcc_R put_R thal_R VP_R amyg_L caud_L GP_L hip_L NAcc_L put_L thal_L VP_L hypothal midbrain brainstem;do
	c=$(($c+1))
	n=`printf %05d $c`
	antsApplyTransforms -d 3 -i ${atldir}/${tiss}.nii.gz -t T1_to_DWI.txt [${sub}_T1rot2atl.txt,1] ${warpdir}/${sub}_rot2atlWarp.nii.gz ${warpdir}/${sub}_rot2atlAffine.txt -r T1_on_DWI.nii.gz -o parcellation/subcorts_full/reg_${n}.nii.gz
done
antsApplyTransforms -d 3 -i ${atldir}/MY_Subcorts_mask.nii.gz -t T1_to_DWI.txt [${sub}_T1rot2atl.txt,1] ${warpdir}/${sub}_rot2atlWarp.nii.gz ${warpdir}/${sub}_rot2atlAffine.txt -r T1_on_DWI.nii.gz -o parcellation/subcorts_full/SubCort_mask.nii.gz
fslmaths parcellation/subcorts_full/SubCort_mask.nii.gz -thr .5 -bin -mas segmentation/GMdeep parcellation/subcorts_full/SubCort_mask.nii.gz

echo "combining subcortical region probabilities into parcellation"
#recombine split files into 4D file
fslmerge -t parcellation/subcorts_full/reg_4D parcellation/subcorts_full/reg_0*nii.gz
#convert 4D file into ranks (and add 1 because fsl uses 0 starting index)
fslmaths parcellation/subcorts_full/reg_4D -rank -Tmaxn -add 1 -mas parcellation/subcorts_full/SubCort_mask.nii.gz parcellation/subcorts_full/SubCort_Parc
#clean up
rm parcellation/subcorts_full/reg_0*nii.gz parcellation/subcorts_full/reg_4D.nii.gz

echo "getting WMGM boundary for subcortical parcellation"
fslmaths parcellation/subcorts_full/SubCort_Parc -mas segmentation/WM_dil parcellation/subcorts_edge/SubCort_Parc

echo "re-split into individual binarized regions for full and edge regions"
c=0
for tiss in Amyg_R caudate_R GP_R HC_R NAcc_R putamen_R thalamus_R VP_R Amyg_L caudate_L GP_L HC_L NAcc_L putamen_L thalamus_L VP_L hypothal midbrain brainstem_ceiling;do
	c=$(($c+1))
	fslmaths parcellation/subcorts_full/SubCort_Parc -thr $c -uthr $c -bin parcellation/subcorts_full/${tiss}
	fslmaths parcellation/subcorts_edge/SubCort_Parc -thr $c -uthr $c -bin parcellation/subcorts_edge/${tiss}
done

;&
prep_probtrack)
mkdir -p probtrackX/network_wholebrain_WMedge
echo -e "\n$sub - prepping probtrackx2"
#echo "SEEDS WILL BE THE WMGM_INNEREDGE - all regions inflated and masked with WM"
#fslmaths parcellation/RMfinal -add parcellation/subcorts_full/SubCort_Parc -bin -dilD -mas segmentation/WM probtrackX/WMGM_inneredge
#echo "SEEDS WILL BE THE WMGM_INNEREDGE - all cort and subcort GM inflated and masked with WM"
#fslmaths segmentation/GMcort -add segmentation/GMdeep -bin -dilD -mas segmentation/WM probtrackX/WMGM_inneredge
echo "SEEDS WILL BE THE WHOLE OF WM"
echo "STOP MASK WILL BE THE EDGE REGIONS"
fslmaths parcellation/RMfinal_edge -add parcellation/subcorts_edge/SubCort_Parc -bin probtrackX/stop_mask
echo "EXCLUSION MASK WILL BE EVERYTHING OTHER THAN EDGE REGIONS OR WM"
fslmaths probtrackX/stop_mask -add segmentation/WM -bin -mul -10 -add segmentation/hardseg_mask -bin probtrackX/exclusion_mask
echo "TARGETS WILL BE THE EDGE REGIONS, order specified in reference file"
cp ${cwd}/probtrack_targetsfile.txt parcellation/targets.txt

#echo "SETUP COMPLETE. RUNNING probtrackx2, seeding from WM-GM inner edge"
#probtrackx2 -x probtrackX/WMGM_inneredge -V 2 -S 600 --steplength=0.3 -P 1000 --forcedir --opd -s bedpostX.bedpostX/merged -m segmentation/hardseg_mask  --dir=probtrackX/network_wholebrain_WMedge --target3=parcellation/targets.txt --avoid=probtrackX/exclusion_mask --omatrix3 --stop=probtrackX/stop_mask

echo "SETUP COMPLETE. RUNNING probtrackx2, seeding from whole WM"
probtrackx2 -x segmentation/WM -V 2 -S 600 --steplength=0.3 -P 2000 --forcedir --opd -s bedpostX.bedpostX/merged -m segmentation/hardseg_mask  --dir=probtrackX/network_wholebrain_WMwhole --target3=parcellation/targets.txt --avoid=probtrackX/exclusion_mask --omatrix3 --stop=probtrackX/stop_mask

esac
popd
exit
