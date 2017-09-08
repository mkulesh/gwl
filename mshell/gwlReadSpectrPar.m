function [aParams]=gwlReadSpectrPar(aFid)
% [aParams]=gwlReadSpectrPar(aFid)
% This procedure reads a binary stream with format 'WP1.4'.
% In GWL do not exist files of 'WP1.4' format, such stream is only a part of
% other GWL binary files, for example, a part of spectrum container 'SP1.3'.
% Therefore, this procedure is only used inside of GWL mshell library.
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%
% Output parameters:   
%   aParams contains the technical details about the used wavelet
%
% This file is part of the GWL library. Copyright (C) 2006-2008 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'WP1.4';
aPar.aObjVer = fread(aFid,5,'*char')';
aPar.aDataType = fread(aFid,1,'uint');
if strcmp(aPar.aObjVer,aVer)==0
    fclose(aFid);
    error(['version or format of input file is incorrect in procedure: gwlReadSpectrPar(). Procedure needs ver. ' aVer]);
end;

aPar.aTransType = fread(aFid,1,'uint');
aPar.aWavelet = fread(aFid,1,'uint');
aPar.aWaveletPar = fread(aFid,1,'*double');
aPar.aIWavelet = fread(aFid,1,'uint');
aPar.aIWaveletPar = fread(aFid,1,'*double');
[aFreq,aParFreq] = gwlReadAxis(aFid);
%aPar.aTMin = fread(aFid,1,'*double');
%aPar.aTMax = fread(aFid,1,'*double');
aNameSize = fread(aFid,1,'uint');
aPar.aWPName = fread(aFid,aNameSize,'*char')';
aPar.aFreq = aFreq;
aPar.aParFreq = aParFreq;
aParams = aPar;

