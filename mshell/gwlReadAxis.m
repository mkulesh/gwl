function [aAxis,aParams]=gwlReadAxis(aFid)
% [aAxis,aParams]=gwlReadAxis(aFid)
% This procedure reads a stream with format 'AX1.4' from a binary data file of GWL.
% The file with this format is produces in the module 'bin/gwlCreateAxis'
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%
% Output parameters:   
%   aAxis is an array contained axis values 
%   aParams contains the technical details about the aAxis variable
%
% Examples: To produce and read an axis, use the following code:
%   gwlExec('gwlCreateAxis',' --count=100 --min=1 --max=1');
%   fid = fopen(aFile,'r');
%   [aAxis,aParams]=gwlReadAxis(fid);
%   fclose(fid);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'AX1.4';
aPar.aObjVer = fread(aFid,5,'*char')';
aPar.aDataType = fread(aFid,1,'uint');
if strcmp(aPar.aObjVer,aVer)==0
    fclose(aFid);
    error(['version or format of input file is incorrect in procedure: gwlReadAxis(). Procedure needs ver. ' aVer]);
end;

aPar.aAxisType = fread(aFid,1,'uint');
aPar.aAxisSign = fread(aFid,1,'uint');
aPar.aMin = fread(aFid,1,'*double');
aPar.aMax = fread(aFid,1,'*double');
aPar.aDelta = fread(aFid,1,'*double');
aPar.aSample = fread(aFid,1,'*double');
[aAxis,aParVec] = gwlReadVector(aFid);

if nargout == 2
     aParams = aPar;
     aParams.aSize = aParVec.aSize;
     aParams.aName = aParVec.aName;
 end;
