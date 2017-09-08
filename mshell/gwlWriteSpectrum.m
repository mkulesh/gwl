function gwlWriteSpectrum(aFid,aTime,aFreq,aSpec,aParams)
% gwlWriteSpectrum(aFid,aTime,aFreq,aSpec,aParams)
% This procedure writes a spectrum as a binary stream with format 'SP1.3'.
% The file with this format is produces also in the module 'bin/gwlCwt'
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%   aTime is an array contained time values
%   aFreq is an array contained frequency values
%   aSpec is a matrix contained spectrum values (real or complex)
%   aParams contains the technical details about the aSpec variable
%
% This file is part of the GWL library. Copyright (C) 2006-2008 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'SP1.3';
aDataType = 8;
if isreal(aSpec) == 0
   aDataType = 16;
end;
fwrite(aFid,aVer,'char')';
fwrite(aFid,aDataType,'uint');

gwlWriteAxis(aFid,aTime,aParams.aParTime);
gwlWriteAxis(aFid,aFreq,aParams.aParFreq);
[aSize1,aSize2,aChanels] = size(aSpec);
fwrite(aFid,aChanels,'uint');
for k=1:aChanels
    gwlWriteMatrix(aFid,aSpec(:,:,k),'spectrum channel');
end;
aNameSize = length(aParams.aName);
fwrite(aFid,aNameSize,'uint');
fwrite(aFid,aParams.aName,'char')';
gwlWriteSpectrPar(aFid,aParams.aParWav);
