%% Convert MATLAB (.FIG) into .png/.jpeg format 
%
%
% f
%
%
% (Nand, 09/2018)

function [] = convertMATLABfig(file, format)

%file  == '...\Brain.fig' 
%file is your MATLAB figure (.fig)

%format == .png, .jpeg, other
%format is the format you want to convert to 

%Do some argument checking
if nargin <1 || nargin > 2 
    error('2 arguments required: file, and format') 
elseif nargin == 1 
    format = '.png';  %Default convert format 
end 

if ~exist(file)
    error('This file does not exist: %s', file)
end 


%Open the MATLAB figure(.fig)
figs = openfig(file);

%Eliminate the .fig handle on MATLAB figure
MAT_fig = file;
eliminate_figextension = find(MAT_fig(:) == '.fig');
MAT_fig(eliminate_figextension(1):end) = []; %eliminate .fig extension

%Add converted format 
convertFIG  = [MAT_fig, format];
saveas(figs(1), convertFIG); %Save the figure into desired format 

if ~exist(convertFIG)
    error('Unable to convert %s into %s format \n', file, format) %Failed Convert 
else
    fprintf('\n%s succesfully converted into %s \n', file, convertFIG)
    Status = 1 %Convert Success!  
end 

end 





