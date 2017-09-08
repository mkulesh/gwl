function [aAxis,aSignal,aParams]=gwlReadSignal(aFid)
% [aAxis,aSignal,aParams]=gwlReadSignal(aFid)
% This procedure reads a signal stream or a binary data file of GWL with format 'SI1.3'.
% The file with this format is produces in modules 'gwlSignalGen', 'gwlSignalRead', 
% 'gwlSignalSum', 'gwlSignalFilter', 'gwlCft', 'gwlIwt', 'gwlAutoCorr', 'gwlConvert'.
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%
% Output parameters:   
%   aAxis is an array contained axis values
%   aSignal is a matrix contained signal values (real or complex)
%   aParams contains the technical details about the aSignal variable
%
% Examples: To read a signal file, use the following code
%   fid = fopen('signal.dat','r'); 
%   [aAxis,aSignal,aParams]=gwlReadSignal(fid); 
%   fclose(fid);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'SI1.3';
aPar.aObjVer = fread(aFid,5,'*char')';
aPar.aDataType = fread(aFid,1,'uint');
if strcmp(aPar.aObjVer,aVer)==0
    fclose(aFid);
    error(['version or format of input file is incorrect in procedure: gwlReadSignal(). Procedure needs ver. ' aVer]);
end;

[aAxis,aParAxis] = gwlReadAxis(aFid);
aPar.aChanCount = fread(aFid,1,'uint');
for k=1:aPar.aChanCount
    [aSignal(:,k),aParVec] = gwlReadVector(aFid);
end;
aNameSize = fread(aFid,1,'uint');
aPar.aName = fread(aFid,aNameSize,'*char')';
aPar.aPointCount = aParAxis.aSize;

if nargout == 3
    aParams = aPar;
    aParams.aMin = aParAxis.aMin;
    aParams.aMax = aParAxis.aMax;
    aParams.aDelta = aParAxis.aDelta;
    aParams.aSample = aParAxis.aSample;
end;
