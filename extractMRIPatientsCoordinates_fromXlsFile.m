% Extracts MRI MNI brain coordinates from patients and areas with "good" cells automatically
%
% cellCounts.counts needs to contain column vector for each session/patient
% requires the following columns in the excel file: 
% Patient ID, MNI areas(left and right Amy, HF, ACC, SMA)
%
% Patient ID and MNI areas have to be type text
% Edited 5/2017 by S. Sullivan
%
%
function [MasterCoords, ALL_MNI2,ALL_MNI_Addl] = extractMRIPatientsCoordinates_fromXlsFile(Excelfile, basepath, cellCountsFile,...
    TargetBrainAreas, BrainAreas, range, AddAnalysis)


%% Section 1: Load info w/ patients and units in brain area for given experiment
%basepath='/home/urut/Dropbox/manuscripts/errorProcessing2015/';
load([basepath, cellCountsFile])
goodcells = {};
CombineSessions = {};
if exist('TEST_cellCounts') ==1
 cellCounts = TEST_cellCounts;
end
for z = 1:length(cellCounts)
    cellCounts_adj(z).sessionName = cellCounts(z).sessionName;
    cellCounts_adj(z).counts = cellCounts(z).counts(:,1);
    if length(cellCounts(z).counts(:,1)) == 8 %no OFC electrodes included
    cellCounts_adj(z).counts(9:10,1) = 0;
    end
end
cellCounts = cellCounts_adj;
% rename patient IDs and keep units for target brain areas
for p = 1:length(cellCounts)
    if ~isempty(strfind(cellCounts(p).sessionName,' ')) == 1
        cellCounts(p).sessionName = strrep(cellCounts(p).sessionName,' ','')
    end
    if isempty(strfind(cellCounts(p).sessionName(1),'P')) == 1 % patient IDs that don't start with P*
        cellCounts(p).sessionName = cellCounts(p).sessionName(1:2);
    elseif ~isempty(strfind(cellCounts(p).sessionName,'CS')) == 1 % CS patients
        nrCharsToUse = 5;
        if length ( cellCounts(p).sessionName ) < 5
            nrCharsToUse = length ( cellCounts(p).sessionName );
        end
        cellCounts(p).sessionName = cellCounts(p).sessionName(1:nrCharsToUse);
    elseif ~isempty(strfind(cellCounts(p).sessionName,'HMH')) == 1 % huntington patients
        nrCharsToUse = 6;
        if length ( cellCounts(p).sessionName ) < 6
            nrCharsToUse = length ( cellCounts(p).sessionName );
        end
        cellCounts(p).sessionName = cellCounts(p).sessionName(1:nrCharsToUse);
    end
    cellCounts(p).counts = cellCounts(p).counts;
    goodcells{1,p} = [cellCounts(p).sessionName];
    goodcells{2,p} = [cellCounts(p).counts];
    
end
for g = 1:length(goodcells)
    if nnz(goodcells{2,g}(TargetBrainAreas)) == 0 % good units found in that area
        goodcells{1,g} = [];
        goodcells{2,g} = [];
    end
end 
goodcells = reshape(goodcells(~cellfun(@isempty,goodcells)),2,[]);

% compiling all sessions for each patient
CombineSessions = unique(goodcells(1,:));
for n = 1:length(CombineSessions)
    for g = 1:length(goodcells)
        if strcmp(goodcells{1,g},CombineSessions{1,n}) == 1
            CombineSessions{2,n}{g} = [goodcells{2,g}] ;
            CombineSessions{2,n} = reshape(CombineSessions{2,n}(~cellfun(@isempty,CombineSessions{2,n})),1,[]);
        end
    end
end
for n = 1:length(CombineSessions)
CombineSessions{2,n} = cell2mat(CombineSessions{2,n});
CombineSessions{2,n} = sum(CombineSessions{2,n},2);
end
clear goodcells
goodcells = CombineSessions;


%% Section 2: Extract Target Coordinates from Patients & Areas with Good Units from Excel Spreadsheet
%basepath='W:\MRIs\';
xlsFile=[basepath Excelfile];
columnPatient=1;% patient ID

% assigning location of Excel Coordinates to brain area based on given experiment 
for b = 1:size(BrainAreas,2)
AllAreas(b) = BrainAreas{2,b};
end
BrainAreas = BrainAreas(1:2,:);
BrainAreas = cell2mat(BrainAreas);

[num,txt,raw] = xlsread(xlsFile, 1, '','basic');  % range selection does not work in basic mode (unix)

MasterCoords=cell(length(goodcells)+1,length(AllAreas)+1); 
MasterCoords{1,1} = 'Patient ID';
AllPt_MasterCoords=cell(length(goodcells)+1,9); 
AllPt_MasterCoords{1,1} = 'Patient ID';

for k = range  %  rows in this worksheet containing all patients
    AllPt_MasterCoords{k-2,1} = raw{k,columnPatient}; % all patients (regardless of "good" units) 
    for a = 1:length(TargetBrainAreas)
         AllPt_MasterCoords{k-2,TargetBrainAreas(a)+1} = raw{k,AllAreas(find(BrainAreas(1,:)==TargetBrainAreas(a)))};
    end
    for g = 1:length(goodcells)
        if ~isempty(strfind(raw{k,columnPatient},goodcells{1,g})) == 1
            MasterCoords{g+1,1} = raw{k,columnPatient};
            for a = 1:length(TargetBrainAreas)   
                if goodcells{2,g}(TargetBrainAreas(a)) > 0 % good units found in that area
                    MasterCoords{g+1,TargetBrainAreas(a)+1} = raw{k,AllAreas(find(BrainAreas(1,:)==TargetBrainAreas(a)))};
                end
            end
        end
    end
end

% excluding coordinates: 
% if no microelectrodes present; if electrode in corpus callosum or not in desired region
for r = 1:size(MasterCoords,1)
    for c = 1:size(MasterCoords,2)
        if isnan(MasterCoords{r,c}) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'callosum')) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'not')) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'maybe')) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'in')) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'out')) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'n/a')) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'no')) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'SEE')) == 1 | ~isempty(strfind(MasterCoords{r,c}, 'too')) == 1 |~isempty(strfind(MasterCoords{r,c}, 'error')) == 1           
            disp(['Removing from coordinate list: ' MasterCoords{r,1} ' because coord is ' num2str(MasterCoords{r,c}) ]);
            MasterCoords{r,c} = []; 
        end
    end
    if size(MasterCoords,2) < length(AllAreas)+1
        MasterCoords(r,length(AllAreas)+1) = [];
    end
end

%all patients & electrodes
for r = 1:size(AllPt_MasterCoords,1)
    for c = 1:size(AllPt_MasterCoords,2)
        if isnan(AllPt_MasterCoords{r,c}) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'callosum')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'not')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'maybe')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'in')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'out')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'n/a')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'no')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'SEE')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'too')) == 1 |~isempty(strfind(AllPt_MasterCoords{r,c}, 'error')) == 1 | ~isempty(strfind(AllPt_MasterCoords{r,c}, 'No')) == 1          
            AllPt_MasterCoords{r,c} = []; 
        end
    end
end

%% Section 3: Organize Target Coordinates
ALL_MNI = {};
for v = 1:length(AllAreas)
    ALL_MNI{1,v} = [MasterCoords(:,v+1)]
    ALL_MNI{1,v} = reshape(ALL_MNI{1,v}(~cellfun(@isempty,ALL_MNI{1,v})),[],1);
    for t = 1:length(ALL_MNI{1,v})
    ALL_MNI2{1,v}(t,1:3) = [str2num(ALL_MNI{1,v}{t})];
    end    
end


%all patients % electrodes - optional
% AllPt_MNI = {};
% for v = 1:length(AllAreas)
%     AllPt_MNI{1,v} = [AllPt_MasterCoords(:,v+1)]
%     AllPt_MNI{1,v} = reshape(AllPt_MNI{1,v}(~cellfun(@isempty,AllPt_MNI{1,v})),[],1);
%     for t = 1:length(AllPt_MNI{1,v})
%     AllPt_MNI2{1,v}(t,1:3) = [str2num(AllPt_MNI{1,v}{t})];
%     end    
% end
% clear ALL_MNI
% ALL_MNI = AllPt_MNI2;

%% Section 4: Additional Analysis
if ~isempty(AddAnalysis) == 1
    if ~isempty(strfind(cellCountsFile, 'stroop')) == 1 % brooks stroop add'l analysis
        goodcellsPts = goodcells(1,:);
        [ALL_MNI_Addl,MasterCoordsAddl, SMA_ErrCongPercent, ACC_ErrCongPercent] = StroopMRI_ErrorCongrCells(AddAnalysis, goodcellsPts, MasterCoords,AllAreas)
    end
else
    ALL_MNI_Addl = {};
end




end