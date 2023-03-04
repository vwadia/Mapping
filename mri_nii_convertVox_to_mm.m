%
% For nii images
% converts from raw voxel coordinates to mm (MNI) coordinates
%
% this code originates in view_nii.m
% 
%urut/march16
function mm = mri_nii_convertVox_to_mm( vox, origin, voxel_size)
mm = ( vox - origin) .* voxel_size;     % this converts the raw coordinates vox to MNI coordinates (mm)
