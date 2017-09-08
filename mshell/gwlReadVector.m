function [aVector,aParams]=gwlReadVector(aFid)
% [aVector,aParams]=gwlReadVector(aFid)
% This procedure reads a binary stream with format 'VC1.3'.
% In GWL do not exist files of 'VC1.3' format, such stream is only a part of
% other GWL binary files, for example, a part of axis container 'AX1.4'.
% Therefore, this procedure is only used inside of GWL mshell library.
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%
% Output parameters:   
%   aVector is an array contained vector values (real or complex)
%   aParams contains the technical details about the aVector variable
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'VC1.3';
aPar.aObjVer = fread(aFid,5,'*char')';
aPar.aDataType = fread(aFid,1,'uint');
if strcmp(aPar.aObjVer,aVer)==0
    fclose(aFid);
    error(['version or format of input file is incorrect in procedure: gwlReadVector(). Procedure needs ver. ' aVer]);
end;

aPar.aSize = fread(aFid,1,'uint');
if aPar.aDataType==8
    aVector(1:aPar.aSize,1) = fread(aFid,aPar.aSize,'*double');
end;
if aPar.aDataType==16
    aCmplVector = fread(aFid,2*aPar.aSize,'*double');
    aVector(1:aPar.aSize,1) = aCmplVector(1:2:(2*aPar.aSize-1)) + i*aCmplVector(2:2:(2*aPar.aSize));
end;
aNameSize = fread(aFid,1,'uint');
aPar.aName = fread(aFid,aNameSize,'*char')';

if nargout == 2
     aParams = aPar;
end;
