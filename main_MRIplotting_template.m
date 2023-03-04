% MRI Plots slices from an MNI brain with markers for electrodes  
%
%==============================PLEASE READ:============================================%
%=================== EDIT SECTION 1 SPECIFIC EXPERIMENT INFORMATION====================%
%=========================THEN, EXECUTE SECTION 2 =====================================%
%========COPY, RENAME, AMEND, & SAVE THIS FILE FOR YOUR SPECIFIC EXPERIMENT============%
%======================================================================================%
% Edited 5/2017
%    OFCs added by S. Sullivan
%
%% SECTION 1: Parameters to Set

% Set basepaths: choose between server or dropbox versions
%== Server version
basepathLocations = 'W:\MRIs\';
%choose template path
basepathMRIs = 'W:\software\amygdalaAtlas\v030216\MNI152_700um\'; % path to 700um template
% or
basepathMRIs = 'W:\software\amygdalaAtlas\v030216\MNI152\'; % path to 1mm template

%== or Dropbox version
basepathMRIs = '~/Dropbox/manuscripts/errorProcessing2015/anatomy/';  
basepathLocations = '~/Dropbox/manuscripts/errorProcessing2015/';


% name of cellCounts mat-file
% cellCounts is a structure array: field 1 = patient/sessions, 
% field 2 = column vector with # of units in each area [0;0;10;3;12;0;9;1]
% for each patient
cellCountsFile =  'cellCountsAllSessions_stroopPaper.mat'; % open this in \\rutishauserulab\Share\MRIs to see example

% OPTIONAL: name of mat-file for additional analysis (# error cells, etc)
% file must be in same format as cellCounts (structure array with 2 fields:
% patients, counts) or another format of your choosing
% this file will be used in customized function that you create later. 
AddAnalysisFile = ' ';



% set which template to use
%   1 = CIT168 MNI 1mm, 2 = CIT168 MNI 700um
templateVersion = 2;


% cellCounts.counts order of brain areas in field 2: 
% e.g. 1=RH, 2=LH, 3=RA, 4=LA, 5=RAC, 6=LAC, 7=RSMA, 8=LSMA, 9=ROFC,
% 10=LOFC
% also assign area and colors for plot --> e.g 'bo' = blue circle, 'co' = cyan circle, etc
RH = 1;   col_RH = 'yo';
LH = 2;   col_LH = 'yo';
RA = 3;   col_RA = 'mo';
LA = 4;   col_LA = 'mo';
RAC = 5;  col_RAC = 'bo'; 
LAC = 6;  col_LAC = 'bo';
RSMA = 7; col_RSMA = 'ro';
LSMA = 8; col_LSMA = 'ro';
ROFC = 9; col_ROFC = 'go'; % may need to change color depending on how well it shows on the figure
LOFC = 10; col_LOFC = 'go';

%Target areas to plot: e.g. [RAC, LAC, ROFC, LOFC]
TargetAreas = [RAC,LAC, RSMA, LSMA, RH, LH, RA, LA];

% Excel File
Excelfilename = 'MNI_electrodePositions.xlsx';
columnPatient=1;% patient ID
ptrange=[6:49]; % indicate rows with patients in column 1, will need to update as more patients are added





%% SECTION 2: EXECUTE THIS WHOLE SECTION TO GENERATE PLOTS (Ctrl + Enter)

% brain area -> excel column
Amy_L = 14;  % left Amygdala coordinates
Amy_R = 23;  % right Amygdala coordinates
HF_L = 34;   % left Hippocampus coordinates
HF_R = 43;   % right Hippocampus coordinates
ACC_L = 52;  % left Anterior Cingulate Cortex coordinates
ACC_R = 62;  % right Anterior Cingulate Cortex coordinates
SMA_L = 71;  % left Supplmentary Motor Area coordinates
SMA_R = 81;  % right Supplementary Motor Area coordinates
OFC_L = 91;    % left Orbitofrontal Cortex coordinates, value may need to change based on Excel File
OFC_R = 100;    % right Orbitofrontal Cortex coordinates

params = [];
params.BrainAreas = {RH,LH,RA,LA,RAC,LAC,RSMA,LSMA,ROFC,LOFC;HF_R,HF_L,Amy_R,Amy_L,ACC_R,ACC_L,SMA_R,SMA_L...
    OFC_R,OFC_L; col_RH,col_LH, col_RA, col_LA, col_RAC, col_LAC, col_RSMA,col_LSMA,col_ROFC,col_LOFC};

params.TargetBrainAreas = TargetAreas;  
params.TemplateVersion = templateVersion;
params.basepathMRIs  = basepathMRIs;
params.basepathLocations = basepathLocations;
params.cellCountsFile = cellCountsFile;
params.Excelfile = Excelfilename;
params.range = ptrange;

if templateVersion == 1
    fnameTemplate = 'CIT168_T1w_1mm_MNI.nii.gz';
elseif templateVersion == 2
    fnameTemplate = 'CIT168_T1w_700um_MNI.nii.gz';
end
params.fnameTemplate = fnameTemplate;

if ~exist('AddAnalysisFile', 'var') || isempty(strfind(AddAnalysisFile, 'mat')) == 1 
    params.AddAnalysis = [];
else
    params.AddAnalysis = AddAnalysisFile;
end

generic_main_plotElectrodePos_onMNI(params)




