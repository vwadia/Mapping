%
% produce cellCounts file for NO data release for MRI plotting
%

clear; clc;

basepaths='V:\LabUsers\carlsona1\Scientific Data Paper\events\';

[NOsessions] = defineNOsessions_SciData();


cellCountsrow = 1
%Load brainArea





for i=1:length(NOsessions)
    if isempty(NOsessions(i).session)
        continue
    end
    brainAreafile = [basepaths NOsessions(i).session '\' NOsessions(i).taskDescr '\brainArea.mat'];
    if ~exist(brainAreafile,'file')
        warning(['didnt find brain area, need to check in dir ' brainAreafile ]);
        continue;
    end
    
    disp(['Processing: ' brainAreafile]);
    
    load(brainAreafile);
    
    %PtID = [NOsessions(i).session];
    
    %Looks at each row in fourth column of brainArea file and determines
    %whether the desired areas are there. If they are, it is noted by a 1
    %in the cellCounts variable where each column is a different area
    
    
    patientID = determinePatientID_NOtask( NOsessions(i) );

    
    cellCounts(cellCountsrow).sessionName =patientID;
    
    
    
    %cellCounts{cellCountsrow,1}=PtID;
    %for k=2:5;
    %cellCounts{cellCountsrow,k}=0;
    %end
    
    cellCountsTmp = [ 0; 0; 0; 0];
    
    %[sum(brainArea(:,4)==4)>0 sum(brainArea(:,4)==3)>0 sum(brainArea(:,4)==2)>0 sum(brainArea(:,4)==4)>1] 
    %^^Julien's suggestion for a more efficient approach 
    for j=1:size(brainArea,1)
        if brainArea(j,4)>4
            continue;
        end
        if brainArea(j,4)==1   %Checks for RH
            cellCountsTmp(1) = cellCountsTmp(1) + 1;  
          
        end
        if brainArea(j,4)==2   %Checks for LH
            cellCountsTmp(2) = cellCountsTmp(2) + 1;  
           
        end
        if brainArea(j,4)==3   %Checks for RA
            cellCountsTmp(3) = cellCountsTmp(3) + 1;  
            
        end
        if brainArea(j,4)==4   %Checks for LA 
            cellCountsTmp(4) = cellCountsTmp(4) + 1;  
          
        end
    end
    
    cellCounts(cellCountsrow).counts = cellCountsTmp;

    cellCountsrow=cellCountsrow+1;
end

cellCountsOrig = cellCounts;

%% collapse cellCounts file so that it contains only one row per patient (and not per session)


allSession_names =[];
for k=1:length(cellCountsOrig)
    allSession_names{k} = cellCountsOrig(k).sessionName;
end

cellCounts_collapsed=[];
cellCounts_collapsed_counter = 0;

allSession_names_unique = unique(allSession_names);
for k=1:length(allSession_names_unique)   
    sessionName_processing = allSession_names_unique(k);
    
    sessionCounts_added=[0; 0; 0; 0];
    for j=1:length(cellCountsOrig)
        if strcmp( cellCountsOrig(j).sessionName, sessionName_processing)==1
            %sessionName_rows = [ sessionName_rows j];
            
            sessionCounts_added = sessionCounts_added + cellCountsOrig(j).counts;
        end
    end
    
    cellCounts_collapsed_counter = cellCounts_collapsed_counter+1;
    
    cellCounts_collapsed(cellCounts_collapsed_counter).sessionName = char(sessionName_processing);
    
    cellCounts_collapsed(cellCounts_collapsed_counter).counts = sessionCounts_added;
end



%%

cellCounts = cellCounts_collapsed;

save('V:\LabUsers\carlsona1\Scientific Data Paper\Scripts\cellCounts_SciDataPaper.mat','cellCounts')


