#!/bin/bash

echo "using the final parcellation MY_Subcorts_BOTH.nii.gz, individual regions are binarized and extracted"

echo "binary mask of all regions combined"
fslmaths MY_Subcorts_BOTH -bin MY_Subcorts_mask

echo "single region regions"
#8 - hypothalamus
#9 - midbrain
#10 - brainstem
fslmaths MY_Subcorts_BOTH.nii.gz -thr 8 -uthr 8 -bin hypothal
fslmaths MY_Subcorts_BOTH.nii.gz -thr 9 -uthr 9 -bin midbrain
fslmaths MY_Subcorts_BOTH.nii.gz -thr 10 -uthr 10 -bin brainstem

echo "two-region regions (uses the corresponding mask accordingly)"
#1 - amygdala
#2 - caudate
#3 - thalamus
#4 - putamen
#5 - hippocampus
#6 - NAcc
#7 - Ventral Pallidum
#11 - GP

for hemi in L R;do
fslmaths MY_Subcorts_BOTH.nii.gz -thr 1 -uthr 1 -bin -mas ../flip_parcs/mask_hardseg${hemi}.nii.gz amyg_${hemi}
fslmaths MY_Subcorts_BOTH.nii.gz -thr 2 -uthr 2 -bin -mas ../flip_parcs/mask_hardseg${hemi}.nii.gz caud_${hemi}
fslmaths MY_Subcorts_BOTH.nii.gz -thr 3 -uthr 3 -bin -mas ../flip_parcs/mask_hardseg${hemi}.nii.gz thal_${hemi}
fslmaths MY_Subcorts_BOTH.nii.gz -thr 4 -uthr 4 -bin -mas ../flip_parcs/mask_hardseg${hemi}.nii.gz put_${hemi}
fslmaths MY_Subcorts_BOTH.nii.gz -thr 5 -uthr 5 -bin -mas ../flip_parcs/mask_hardseg${hemi}.nii.gz hip_${hemi}
fslmaths MY_Subcorts_BOTH.nii.gz -thr 6 -uthr 6 -bin -mas ../flip_parcs/mask_hardseg${hemi}.nii.gz NAcc_${hemi}
fslmaths MY_Subcorts_BOTH.nii.gz -thr 7 -uthr 7 -bin -mas ../flip_parcs/mask_hardseg${hemi}.nii.gz VP_${hemi}
fslmaths MY_Subcorts_BOTH.nii.gz -thr 11 -uthr 11 -bin -mas ../flip_parcs/mask_hardseg${hemi}.nii.gz GP_${hemi}
done
