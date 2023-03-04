%
%reads DICOM files that come from the clinical HMH scanner, and writes them
%again, using a compression syntax that is readable by brainvoyager.
%
%urut/sept05
function convertDICOM( origName, outName)

infoOrig = dicominfo(origName);
Y = dicomread(origName);


info.SOPClassUID='1.2.840.10008.5.1.4.1.1.4';

%---- copy only wanted fields; rest gets discarded
fieldsToCopy={
'MRAcquisitionType' 'SequenceName' 'SliceThickness' 'RepetitionTime' 'EchoTime' ...
'InversionTime' 'NumberOfAverages' 'ImagingFrequency' 'ImagedNucleus' 'EchoNumbers' ...
'SpacingBetweenSlices' 'NumberOfPhaseEncodingSteps' 'EchoTrainLength' 'PercentPhaseFieldOfView'...
'DeviceSerialNumber' 'SoftwareVersions' 'ProtocolName' 'AcquisitionMatrix' 'PhaseEncodingDirection'...
'FlipAngle' 'PatientPosition' 'StudyInstanceUID' 'SeriesInstanceUID' 'StudyID' 'SeriesNumber'...
'AcquisitionNumber' 'InstanceNumber' 'PatientOrientation' 'ImagePositionPatient' 'ImageOrientationPatient'...
'FrameOfReferenceUID' 'Laterality' 'NumberOfTemporalPositions' 'PositionReferenceIndicator' 'SliceLocation' 'SamplesPerPixel'...
'ImageComments' 'ProtocolName' 'SoftwareVersions' 'PatientsSex' 'PatientsAge' 'PatientsWeight' 'PatientComments','PatientsName',...
'ScanningSequence' 'PatientID' 'PatientsBirthDate' 'SeriesDescription' 'Modality' 'Manufacturer' 'InstitutionName'...
'StudyDate' 'SeriesDate' 'AcquisitionDate' 'StudyTime' 'SeriesTime' 'AcquisitionTime' 'PixelSpacing'
};

for i=1:length( fieldsToCopy )
    if isfield( infoOrig, fieldsToCopy{i} )
        eval( [ 'info.' fieldsToCopy{i} '= infoOrig.' fieldsToCopy{i} ';'] );    
    end    
end

dicomwrite(Y,outName,info,'TransferSyntax','1.2.840.10008.1.2.1');
