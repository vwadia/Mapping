
basepath='F:\brainvoyager\new\DG\';
basepathOrig = [basepath 'orig\'];
basepathConv = [basepath 'conv\'];
basename = 'view';

medconpath = 'C:\Program Files\XMedCon\bin\medcon.exe';

%HM
%sagittal=1:19;
%horizontal=20:42;
%horizontal2=43:65;
%coronal=66:85;

%RA
sagittal=1:19;
horizontal=20:42;
horizontal2=43:65;
coronal=66:88;

imgsToLoad = [ sagittal horizontal horizontal2 coronal ];
basenameFile=[basepath basename];

for i=1:length(imgsToLoad)
    origName = [basepathOrig '\view' num2str(imgsToLoad(i),'%.4d') '.dcm' ];
    outName = [basepath 'conv\view' num2str(imgsToLoad(i),'%.4d') '.dcm'];
    
    disp(['converting ' origName]);
    convertDICOM( origName, outName);    
end

