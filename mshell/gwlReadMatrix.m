function [aMatrix,aParams]=gwlReadMatrix(aFid)
% [aMatrix,aParams]=gwlReadMatrix(aFid) 
% This procedure reads a binary stream with format 'MT1.3'.
% In GWL do not exist files of 'MT1.3' format, such stream is only a part of
% other GWL binary files, for example, a part of spectrum container 'SP1.3'.
% Therefore, this procedure is only used inside of GWL mshell library.
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%
% Output parameters:   
%   aMatrix is an array contained matrix values (real or complex)
%   aParams contains the technical details about the aMatrix variable
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'MT1.3';
aPar.aObjVer = fread(aFid,5,'*char')';
aPar.aDataType = fread(aFid,1,'uint');
if strcmp(aPar.aObjVer,aVer)==0
    fclose(aFid);
    error(['version or format of input file is incorrect in procedure: gwlReadMatrix(). Procedure needs ver. ' aVer]);
end;

aPar.aRows = fread(aFid,1,'uint');
aPar.aCols = fread(aFid,1,'uint');
aSize = aPar.aRows*aPar.aCols;
if aPar.aDataType==8
    aVector(1:aSize,1) = fread(aFid,aSize,'*double');
end;
if aPar.aDataType==16
    aCmplVector = fread(aFid,2*aSize,'*double');
    aVector(1:aSize,1) = aCmplVector(1:2:(2*aSize-1)) + i*aCmplVector(2:2:(2*aSize));
end;
aMatrix = transpose(reshape(aVector,aPar.aCols,aPar.aRows));
aNameSize = fread(aFid,1,'uint');
aPar.aName = fread(aFid,aNameSize,'*char')';

if nargout == 2
     aParams = aPar;
end;
