function [aTime,aFreq,aSpec,aParams]=gwlReadSpectrum(aFid)
% [aTime,aFreq,aSpec,aParams]=gwlReadSpectrum(aFid)
% This procedure reads a spectrum stream or a binary data file of GWL with format 'SP1.3''.
% The file with this format is produces in modules 'gwlCwt', 'gwlET2D', 
% 'gwlET3D', 'gwlDiffeoLin', 'gwlDiffeoDisp', 'gwlTransFK'.
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%
% Output parameters:   
%   aTime is an array contained time values
%   aFreq is an array contained frequency values
%   aSpec is a matrix contained spectrum values (real or complex)
%   aParams contains the technical details about the aSpec variable
%
% Examples: To read a spectrum file, use the following code
%   fid = fopen('spectrum.dat','r'); 
%   [aAxis,aSignal,aParams]=gwlReadSpectrum(fid); 
%   fclose(fid);
%
% This file is part of the GWL library. Copyright (C) 2006-2008 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'SP1.3';
aPar.aObjVer = fread(aFid,5,'*char')';
aPar.aDataType = fread(aFid,1,'uint');
if strcmp(aPar.aObjVer,aVer)==0
    fclose(aFid);
    error(['version or format of input file is incorrect in procedure: gwlReadSpectrum(). Procedure needs ver. ' aVer]);
end;

[aTime,aParTime] = gwlReadAxis(aFid);
[aFreq,aParFreq] = gwlReadAxis(aFid);
aPar.aChanels = fread(aFid,1,'uint');
for k=1:aPar.aChanels
    [aSpec(:,:,k),aParMat] = gwlReadMatrix(aFid);
end;
aNameSize = fread(aFid,1,'uint');
aPar.aName = fread(aFid,aNameSize,'*char')';
if nargout == 4
     aParams = aPar;
     aParams.aVoices = aParMat.aRows;
     aParams.aPoints = aParMat.aCols;
     aParams.aParTime = aParTime;
     aParams.aParFreq = aParFreq;
     aParams.aParWav = gwlReadSpectrPar(aFid);
end;

