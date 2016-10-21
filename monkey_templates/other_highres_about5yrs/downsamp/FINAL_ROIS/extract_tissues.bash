#!/bin/bash

echo "using the final segmentation HardSeg_BOTH.nii.gz, individual tissue clases are binarized and extracted"

fslmaths HardSeg_BOTH.nii.gz -thr 1 -uthr 1 -bin tissue_CSF
fslmaths HardSeg_BOTH.nii.gz -thr 2 -uthr 2 -bin tissue_GMcort
fslmaths HardSeg_BOTH.nii.gz -thr 3 -uthr 3 -bin tissue_WM
fslmaths HardSeg_BOTH.nii.gz -thr 4 -uthr 4 -bin tissue_GMdeep

fslmaths HardSeg_BOTH.nii.gz -bin HardSeg_mask.nii.gz

exit
