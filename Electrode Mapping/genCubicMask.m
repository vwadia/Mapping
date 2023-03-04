% By Shuo Wang Jan 10 2013
% generate cubic mask

clear; clc;

patientID = 'test';

dataInDir = '/Volumes/adolphsusers/ueli/MRI/Masks1/';
dataOutDir = '/Volumes/adolphsusers/ueli/MRI/Masks/';

MNI = [182 218 182]; % x, y, z
offset = [90 -126 -72];
r = 2;

switch patientID
    case 'p17'
        xL = -22; yL = -2; zL = -24;
        xR = 18; yR = 0; zR = -22;
    case 'p21'
        xL = -20; yL = 0; zL = -20;
        xR = 16; yR = -2; zR = -18;
    case 'p25'
        xL = -21; yL = -5; zL = -25;
        xR = 18; yR = -5; zR = -24;
    case 'p27'
        xL = -15; yL = -2; zL = -22;
        xR = 18; yR = -2; zR = -22;
    case 'p28'
        xL = -18; yL = -4; zL = -24;
        xR = 20; yR = -4; zR = -24;
    case 'p29'
        xL = -18; yL = -6; zL = -18;
        xR = 20; yR = -2; zR = -22;
    case 'test'
        xL = -22; yL = -2; zL = -24;
        xR = 18; yR = 0; zR = -22;
        c = 20;
end

xL = offset(1) - xL + 1;
yL = yL - offset(2) + 1;
zL = zL - offset(3) + 1;

xR = offset(1) - xR + 1;
yR = yR - offset(2) + 1;
zR = zR - offset(3) + 1;

%% load Nii

if ~isempty(dir([dataInDir filesep patientID '.nii']))
    Nii = load_nii([dataInDir filesep patientID '.nii']);
    m = Nii.hdr.dime.glmax;
else
    dataInDir = '/Volumes/adolphsusers/ueli/MRI/Mask Templates/';
    Nii = load_nii([dataInDir filesep 'MNI152_T1_1mm-mask.nii']);
    m = c;
end


%% generate mask

M = zeros(MNI);

for i = 1:MNI(1)
    for j = 1:MNI(2)
        for k = 1:MNI(3)
            if i>xL-r && i<xL+r && j>yL-r && j<yL+r && k>zL-r && k<zL+r
                M(i,j,k) = m;
            end
        end
    end
end

disp(sum(M(:)>0));

for i = 1:MNI(1)
    for j = 1:MNI(2)
        for k = 1:MNI(3)
            if i>xR-r && i<xR+r && j>yR-r && j<yR+r && k>zR-r && k<zR+r
                M(i,j,k) = m;
            end
        end
    end
end

disp(sum(M(:)>0));

%% save output

M = uint16(M);

Nii.img = M;
Nii.hdr.hist.aux_file = 'Random-Rainbow';

save_nii(Nii,[dataOutDir filesep patientID '.nii']);
