%
%returns transformation matrix for a dicom file
%
%urut/sept05
function [M] = getHeaderInfo( origName )

info = dicominfo(origName);

%notation according to 7.6.2.1.1 of dicom 3.0 standard

%0020,0032
imgPos = info.ImagePositionPatient;   
Sxyz = imgPos;

%0020,0037
imgOrient = info.ImageOrientationPatient;  
Xxyz = imgOrient(1:3);
Yxyz = imgOrient(4:6);

%0028,0030. see table C.7-10
spacing = info.PixelSpacing;
dj = spacing(1); %row spacing
di = spacing(2); %column spacing

%transformation into RCS space
M = [Xxyz(1)*di Yxyz(1)*dj 0 Sxyz(1);
     Xxyz(2)*di Yxyz(2)*dj 0 Sxyz(2);
     Xxyz(3)*di Yxyz(3)*dj 0 Sxyz(3);
         0          0      0    1;
    ];
    
