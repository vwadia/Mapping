%% plot slices from an MNI brain with dots superimposed
%
% version that can use either 700um or 1mm template
%
%

function ALL_MNI2 = generic_main_plotElectrodePos_onMNI(params)

%% Section 1: load MNI template brain
% if you want a different slice, change the variables (sliceNrX/Y/Z) below
switch( params.TemplateVersion )
    case 1  %CIT168 MNI 1mm
        
        sliceNrX = 100;  %  x - medial to lateral.  A saggital plane
        sliceNrY = 121;    % y - anterior to posterior. A corronal plane
        sliceNrZ = 80;    % z - ventral to dorsal (up/down). An axianl plane
        
        xlims_sag = [20 200];
        ylims_sag = [10 150];
        
    case 2  %CIT168 MNI 700um
        
%Hints in Plotting: 
%SMA: Saggital (use x = 122), Axial (use z = 165), Coronal (use y = 209)
%ACC: Saggital (use x = 125), Axial (use z = 138 ), Coronal (y = 209)
%OFC: Saggital (use x = 125 ), Axial (use z = 90 ), Coronal (y =  )
%Amy: Saggital (use x = 163), Axial (use z = 81), Coronal (y = 177)
%Hip: Saggital (use x = ), Axial (use z =), Coronal (y = )

        sliceNrX = 122;  %  x - medial to lateral.  A saggital plane  (for hipp/amy elec, try 163)
        sliceNrY = 209;    % y - anterior to posterior. A corronal plane
        sliceNrZ = 165;    % z - ventral to dorsal (up/down). An axial plane (for hipp/amy elec, try 82)
         
        xlims_sag = [40 280];
        ylims_sag = [10 220];
%         xlims_sag = [180 280]; % starts at origin, only ACC/SMA region shown
%         ylims_sag = [100 220]; % starts at origin, only ACC/SMA region shown
end

% which MNI coordinates should be labeled (tickmarks and text)
% ax_forAxial_requested = [-90:10:60];
% ay_forAxial_requested = [-60:10:60];
% 
% ax_forCor_requested = [-40:10:60];
% ay_forCor_requested = [-60:10:60];
% 
% ax_forSag_requested = [-80:10:80];
% ay_forSag_requested = [-60:10:60];

% another option for tickmarks (increments of 20)
ax_forAxial_requested = [-100:20:80];
ay_forAxial_requested = [-80:20:80];

ax_forCor_requested = [-60:20:80];
ay_forCor_requested = [-80:20:80];

ax_forSag_requested = [-100:20:80];
ay_forSag_requested = [-80:20:80];

%% load the template and view it
fname = [params.basepathMRIs params.fnameTemplate];
nii = load_nii(fname);

% which slice numbers are at the origin?
origin = abs(nii.hdr.hist.originator(1:3));   % these are slice numbers. Using these for sliceNrX/Y/Z below will show the slices at MNI (0,0,0)
voxel_size = abs(nii.hdr.dime.pixdim(2:4));		% vol in mm

% convert the desired label coordinates (in MNI) into voxel coordinates
ax_forSag = mri_nii_convertmm_to_vox( ax_forSag_requested,  origin(2), voxel_size(2) );
ay_forSag = mri_nii_convertmm_to_vox( ay_forSag_requested,  origin(3), voxel_size(2) );

ax_forAxial = mri_nii_convertmm_to_vox( ax_forAxial_requested,  origin(2), voxel_size(1) );
ay_forAxial = mri_nii_convertmm_to_vox( ay_forAxial_requested,  origin(1), voxel_size(1) );

ax_forCor = mri_nii_convertmm_to_vox( ax_forCor_requested,  origin(3), voxel_size(3) );
ay_forCor = mri_nii_convertmm_to_vox( ay_forCor_requested,  origin(1), voxel_size(3) );

%% Section 2: inspect visual (GUI)
option.setunit='mm';
view_nii(nii, option);

%% Section 3: pick a slice and prepare for display
sliceToPlot_Z = nii.img(:,:,sliceNrZ);   %axial
sliceToPlot_Y = squeeze(nii.img(:,sliceNrY,:));   %cor
sliceToPlot_X = squeeze( nii.img(sliceNrX,:,:) );   %sage

%=== what is the MNI coordinate of this frame?  Verify value of mm
sag=sliceNrX;
cor=sliceNrY;
axi=sliceNrZ;
vox = [sag,cor,axi];
mm = mri_nii_convertVox_to_mm( vox, origin, voxel_size)

%% Section 4: plot the selected slices
%axial
figure(26);
colormap('gray');

imagesc(sliceToPlot_Z )   % display a slice
set(gca,'xtick',ax_forAxial,'XColor', [0.6 0.2 0]); 
set(gca,'ytick',ay_forAxial, 'YColor', [0.6 0.2 0]); 
set(gca,'XTickLabel',round(mri_nii_convertVox_to_mm([ax_forAxial], origin(2), voxel_size(1) )));
set(gca,'YTickLabel',round(mri_nii_convertVox_to_mm([ay_forAxial], origin(1), voxel_size(1) )));
title(['Axial z#=' num2str(sliceNrZ) ' MNI z=' num2str(mm(3))], 'FontSize' ,18);
xlabel('MNI y Coordinate (mm)', 'FontSize', 18);
ylabel('MNI x Coordinate (mm)', 'FontSize', 18);
axis equal
axis tight
%coronal
figure(27); 
colormap('gray');
imagesc(sliceToPlot_Y )   % display a slice
set(gca,'xtick',ax_forCor,'XColor', [0.6 0.2 0]);
set(gca,'ytick', ay_forCor, 'YColor', [0.6 0.2 0]); 

set(gca,'XTickLabel',round(mri_nii_convertVox_to_mm([ax_forCor], origin(3), voxel_size(1) )));
set(gca,'YTickLabel',round(mri_nii_convertVox_to_mm([ay_forCor], origin(1), voxel_size(1) )));
title(['Cor y#=' num2str(sliceNrY) ' MNI y=' num2str(mm(2))], 'FontSize' ,18);
xlabel('MNI z Coordinate (mm)', 'FontSize', 18);
ylabel('MNI x Coordinate (mm)', 'FontSize', 18);
axis equal
axis tight

%sagittal
figure(28);
% subplot(2,1,2);
colormap('gray');
sliceToPlot_X_rotated = permute( sliceToPlot_X, [2 1 3]);
imagesc(sliceToPlot_X_rotated )   % display a slice

set(gca,'xtick',ax_forSag,'XColor', [0.6 0.2 0]); 
set(gca,'ytick',ay_forSag, 'YColor', [0.6 0.2 0]); 

set(gca,'XTickLabel',round(mri_nii_convertVox_to_mm([ax_forSag], origin(2), voxel_size(1) )));
set(gca,'YTickLabel',round(mri_nii_convertVox_to_mm([ay_forSag], origin(3), voxel_size(1) )));
axis equal

title(['Sag x#=' num2str(sliceNrX) ' MNI x=' num2str(mm(1))], 'FontSize' ,18);
set(gca,'XDir','normal')
set(gca,'YDir','normal')
xlabel('MNI y Coordinate (mm)', 'FontSize', 18);
ylabel('MNI z Coordinate (mm)', 'FontSize', 18);

if ~isempty(xlims_sag)
    xlim(xlims_sag);
end
if ~isempty(ylims_sag)
    ylim(ylims_sag);
end
%% Section 5: import MNI coordinates and convert to voxel coordinates
% export coordinates from Excel spreadsheet


[MasterCoords, ALL_MNI2,ALL_MNI_Addl] = extractMRIPatientsCoordinates_fromXlsFile(params.Excelfile, params.basepathLocations, params.cellCountsFile,...
    params.TargetBrainAreas, params.BrainAreas, params.range, params.AddAnalysis)


% according to the mean position, pick the best slices to collapse to for the different dimensions
for i = 1:length(ALL_MNI2)
    if ~isempty(ALL_MNI2{i}) == 1
        MeanPos{i} = mean(ALL_MNI2{i});
    end
end

% convert from MNI to voxel coordinates
for i = 1:length(ALL_MNI2)
    if ~isempty(ALL_MNI2{i}) == 1
        MNI_converted{i} = mri_nii_convertmm_to_vox( ALL_MNI2{i}, origin, voxel_size);
    end
end


%% Section 6a): plot MNI coordinates on top of saggital (left brain side)
% subplot(2,2,3);
% edgecolor = {'g','k','g','k','g','k','g','k'}; % right = green, left  = black edge
figure(28)
% ignoring the x axis, plot all dots in y/z plane on saggital slice
% collapse right onto left
hold on
BrainAreas = params.BrainAreas(1,:);
BrainAreas = cell2mat(BrainAreas);
for R = 1:length(MNI_converted)
    if ~isempty(MNI_converted{R}) == 1
        plot( MNI_converted{R}(:,2), MNI_converted{R}(:,3), params.BrainAreas{3,find(BrainAreas ==R)}, 'MarkerSize', 4, 'MarkerFaceColor', params.BrainAreas{3,find(BrainAreas ==R)}(1))%, 'MarkerEdgeColor', edgecolor{R})
    end
end

hold off

%% Section 6b): plot MNI coordinates on top of coronal 
% subplot(2,2,2);
figure(27) 
% plot only those within range (few millimeters)
tolerance = 10;   % +- so many mm
range = [sliceNrY-tolerance sliceNrY+tolerance];  
for R = 1:length(MNI_converted)
    if ~isempty(MNI_converted{R}) == 1
        indsToPlot{R} = find( MNI_converted{R}(:,2)>range(1) & MNI_converted{R}(:,2)<range(2) );
    end
end

hold on
for R = 1:length(MNI_converted)
    if ~isempty(MNI_converted{R}) == 1
        plot( MNI_converted{R}(indsToPlot{R},3), MNI_converted{R}(indsToPlot{R},1), params.BrainAreas{3,find(BrainAreas ==R)}, 'MarkerSize', 4, 'MarkerFaceColor', params.BrainAreas{3,find(BrainAreas ==R)}(1))%, 'MarkerEdgeColor', edgecolor{R})
    end
end
hold off

%% Section 6c): plot MNI coordinates on top of axial/horizontal
% subplot(2,2,1);
figure(26)
tolerance = 10;   % +- so many mm, can edit tolerance if you want more coordinate to be included
range = [sliceNrZ-tolerance sliceNrZ+tolerance];  
for R = 1:length(MNI_converted)
    if ~isempty(MNI_converted{R}) == 1
        indsToPlot{R} = find( MNI_converted{R}(:,3)>range(1) & MNI_converted{R}(:,3)<range(2) );
    end
end

hold on
for R = 1:length(MNI_converted)
    if ~isempty(MNI_converted{R}) == 1
        plot( MNI_converted{R}(indsToPlot{R},2), MNI_converted{R}(indsToPlot{R},1), params.BrainAreas{3,find(BrainAreas ==R)}, 'MarkerSize', 4, 'MarkerFaceColor', params.BrainAreas{3,find(BrainAreas ==R)}(1))%, 'MarkerEdgeColor', edgecolor{R})
    end
end
hold off


end