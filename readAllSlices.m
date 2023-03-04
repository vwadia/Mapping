function [Pall, inds] = readAllSlices( imgsToLoad, basenameFile )

Pall=[];
inds=[]; %from/to that corresponds to a certain img
for kk=1:length(imgsToLoad)
    imgInd = imgsToLoad(kk);
    origName1=[ basenameFile num2str(imgInd,'%.4d') '.dcm'];
    M1 = getHeaderInfo( origName1 );
    img1 = dicomread(origName1);
    P1 = transformImg(img1,M1);
    
    before=size(Pall,1);
    
    Pall = single([Pall; P1]);
    
    inds(kk,1:2) =  [before+1 size(Pall,1)];
    
    disp(['slice loaded: ' origName1]);
end