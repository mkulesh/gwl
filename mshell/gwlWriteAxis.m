function gwlWriteAxis(aFid,aAxis,aParams)
% gwlWriteAxis(aFid,aAxis,aParams)
% This procedure writes an axis as a binary stream with format 'AX1.4'.
% The file with this format is produces also in the module 'bin/gwlCreateAxis'
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%   aAxis is an array contained axis values 
%   aParams contains the technical details about the aAxis variable, it is an
%   output of gwlReadAxis procedure.
%
% This file is part of the GWL library. Copyright (C) 2006-2008 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

fwrite(aFid,aParams.aObjVer,'char')';
fwrite(aFid,aParams.aDataType,'uint');
fwrite(aFid,aParams.aAxisType,'uint');
fwrite(aFid,aParams.aAxisSign,'uint');
fwrite(aFid,aParams.aMin,'double');
fwrite(aFid,aParams.aMax,'double');
fwrite(aFid,aParams.aDelta,'double');
fwrite(aFid,aParams.aSample,'double');
gwlWriteVector(aFid,aAxis,aParams.aName);
