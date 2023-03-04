
% plots figure of structural MRIs


%basepath = 'C:\data\mris\P11\data\EXP00000\';
%patientID='P11';
%imgsToPlot = [39 86 92] +1 ;  %zero-based, as jpg!
%IDs        ={'AX FLAIR','SAG T1 FLAIR left hippo','SAG T1 FLAIR right amy'};

basepath = 'W:\MRIs\mris\P17\data\EXP00000';
patientID='P17';

imgsToPlot = [59 60 61]  ;  %zero-based, as jpg!
IDs        ={'T2','T2','T2'};

for k=1:length(imgsToPlot)
    fname=[basepath '\view' num2str(imgsToPlot(k),'%.4d') '.dcm'];
    
   % infoOrig = dicominfo(fname); %dont need this,only if wanting to look at header
    
    img = dicomread(fname);

    figure(k);
    
    imshow( img );
    imcontrast;
    
    title([ patientID ' ' IDs{k}]);
    
   
    
    
end



