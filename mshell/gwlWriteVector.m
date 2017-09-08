function gwlWriteVector(aFid,aVector,aName)
% gwlWriteVector(aFid,aVector,aName)
% This procedure writes a binary stream with format 'VC1.3'.
% In GWL do not exist files of 'VC1.3' format, such stream is only a part of
% other GWL binary files, for example, a part of axis container 'AX1.4'.
% Therefore, this procedure is only used inside of GWL mshell library.
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'wb') command.
%   aVector is an array contained vector values (only real)
%   aName is the name of the vector
%
% This file is part of the GWL library. Copyright (C) 2006-2008 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'VC1.3';
aDataType = 8;
fwrite(aFid,aVer,'char')';
fwrite(aFid,aDataType,'uint');

aDim = 1;
%fwrite(aFid,aDim,'uint');
aSize = length(aVector);
fwrite(aFid,aSize,'uint');
fwrite(aFid,aVector,'double');
aNameSize = length(aName);
fwrite(aFid,aNameSize,'uint');
fwrite(aFid,aName,'char')';
