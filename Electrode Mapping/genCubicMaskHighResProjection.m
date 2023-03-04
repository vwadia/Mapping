% By Shuo Wang Jan 10 2013
% generate cubic mask
% modified by Shuo Wang on Jan 13 2013 for 0.5 mm resolution
% project along y-axis

clear; clc;

patientID = 'p17';

dataInDir = '/Volumes/adolphsusers/ueli/MRI/Mask Templates/';
dataOutDir = '/Volumes/adolphsusers/ueli/MRI/Masks High Resolution Projection/';

MNI = [182 218 182]*2; % x, y, z
offset = [90 -126 -72];
r = 2;

switch patientID
    case 'p17'
        xL = -22; yL = -2; zL = -24;
        xR = 18; yR = 0; zR = -22;
        c = 1;
    case 'p21'
        xL = -20; yL = 0; zL = -20;
        xR = 16; yR = -2; zR = -18;
        c = 2;
    case 'p25'
        xL = -21; yL = -5; zL = -25;
        xR = 18; yR = -5; zR = -24;
        c = 3;
    case 'p27'
        xL = -15; yL = -2; zL = -22;
        xR = 18; yR = -2; zR = -22;
        c = 4;
    case 'p28'
        xL = -18; yL = -4; zL = -24;
        xR = 20; yR = -4; zR = -24;
        c = 5;
    case 'p29'
        xL = -18; yL = -6; zL = -18;
        xR = 20; yR = -2; zR = -22;
        c = 6;
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

xL = xL * 2;
yL = yL * 2;
zL = zL * 2;

xR = xR * 2;
yR = yR * 2;
zR = zR * 2;

%% load Nii

Nii = load_nii([dataInDir filesep 'MNI152_T1_0.5mm-mask.nii']);
m = c;

%% generate mask

M = zeros(MNI);

for i = 1:MNI(1)
    for j = 1:MNI(2)
        for k = 1:MNI(3)
            if i>xL-r && i<xL+r && k>zL-r && k<zL+r
                M(i,j,k) = m;
            end
        end
    end
end

disp(sum(M(:)>0));

for i = 1:MNI(1)
    for j = 1:MNI(2)
        for k = 1:MNI(3)
            if i>xR-r && i<xR+r && k>zR-r && k<zR+r
                M(i,j,k) = m;
            end
        end
    end
end

disp(sum(M(:)>0));

%% save output

M = int16(M);

Nii.img = M;
Nii.hdr.hist.aux_file = 'Random-Rainbow';

save_nii(Nii,[dataOutDir filesep patientID '.nii']);
