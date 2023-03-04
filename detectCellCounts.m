%% This function will get the cell counts within BrainArea.mat 
%
% Note: The labeling of your cells (in brainArea.mat) must match the format 
% below 
% 1=RH, 2=LH, 3=RA, 4=LA, 5=RAC, 6=LAC, 7=RSMA, 8=LSMA, 11=ROFC, 12=LOFC
% 
% (Nand 09/2018)



function [getCellCounts] = detectCellCounts(file)

%file is the location of "...\brainArea.mat" 

if ~exist(file)
    error('This file does not exist: %s', file) 
else 
    load(file);
end 


%Delete the channels with no units 
toDelete = brainArea(:, 2) < 1;
brainArea(toDelete, :) = [] ;


%Count the number of units within each channel 
% 1=RH, 2=LH, 3=RA, 4=LA, 5=RAC, 6=LAC, 7=RSMA, 8=LSMA, 11=ROFC, 12=LOFC

%Get the number of units in RH (1)         
RH = find(brainArea(:, 4) == 1 );
numUnits_RH = sum(brainArea(RH, 2));
      
%Get the number of units in LH (2)         
LH = find(brainArea(:, 4) == 2 );
numUnits_LH = sum(brainArea(LH, 2));
                      
%Get the number of units in RA (3)         
RA = find(brainArea(:, 4) == 3 );
numUnits_RA = sum(brainArea(RA, 2));

%Get the number of units in LA (4)         
LA = find(brainArea(:, 4) == 4 );
numUnits_LA = sum(brainArea(LA, 2));
                
%Get the number of units in RAC (5)
RAC = find(brainArea(:, 4) == 5 );
numUnits_RAC = sum(brainArea(RAC, 2));

%Get the number of units in LAC (6)
LAC = find(brainArea(:, 4) == 6 );
numUnits_LAC = sum(brainArea(LAC, 2));

%Get the number of units in RSMA (7)
RSMA = find(brainArea(:, 4) == 7 );
numUnits_RSMA = sum(brainArea(RSMA, 2));

%Get the number of units in LSMA (8)
LSMA = find(brainArea(:, 4) == 8 );
numUnits_LSMA = sum(brainArea(LSMA, 2));

%Get the number of units in ROFC (11)
ROFC = find(brainArea(:, 4) == 11 );
numUnits_ROFC = sum(brainArea(ROFC, 2));

%Get the number of units in LOFC (12)
LOFC = find(brainArea(:, 4) == 12 );
numUnits_LOFC = sum(brainArea(LOFC, 2));


total_numofUnits = [numUnits_RH; numUnits_LH; numUnits_RA; numUnits_LA; numUnits_RAC; numUnits_LAC; numUnits_RSMA; numUnits_LSMA; numUnits_ROFC; numUnits_LOFC]; 
    

getCellCounts.brainArea = total_numofUnits;



%Do some validation of cellCounts 
if ~(sum(getCellCounts.brainArea) == sum(brainArea(:, 2)))
    error('Incorrect Cell Counts, Need to Fix Cell Counts')
end 



end 
