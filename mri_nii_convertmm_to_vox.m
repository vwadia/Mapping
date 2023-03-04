%
% For nii images
% converts from mm (MNI) coordinates to raw voxel coordinates
%
% this code originates in view_nii.m
% 
% if mm has more than one row, is applied independenly for each row
%
%
%urut/march16
function vox = mri_nii_convertmm_to_vox( mm, origin, voxel_size)
   
%convert MNI to pixels
for k=1:size(mm,1)
    vox(k,:) = mm(k,:)./voxel_size + origin;
end
