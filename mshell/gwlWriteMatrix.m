function gwlWriteMatrix(aFid,aMatrix,aName)
% gwlWriteMatrix(aFid,aMatrix,aName)
% This procedure writes a binary stream with format 'MT1.3'.
% In GWL do not exist files of 'MT1.3' format, such stream is only a part of
% other GWL binary files, for example, a part of spectrum container.
% Therefore, this procedure is only used inside of GWL mshell library.
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'wb') command.
%   aVector is an array contained matrix values (real or complex)
%   aName is the name of the matrix
%
% This file is part of the GWL library. Copyright (C) 2006-2008 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'MT1.3';
aDataType = 8;
if isreal(aMatrix) == 0
   aDataType = 16;
end;
fwrite(aFid,aVer,'char')';
fwrite(aFid,aDataType,'uint');

aDim = 2;
%fwrite(aFid,aDim,'uint');
[aSize1,aSize2] = size(aMatrix);
fwrite(aFid,aSize1,'uint');
fwrite(aFid,aSize2,'uint');
aSize = aSize1*aSize2;
aVector = reshape(transpose(aMatrix),aSize,1);
if aDataType==8
    fwrite(aFid,aVector,'double');
end;
if aDataType==16
    aCmplVector(1:2:(2*aSize-1)) = real(aVector);
    aCmplVector(2:2:(2*aSize)) = imag(aVector);
    fwrite(aFid,aCmplVector,'double');
end;
aNameSize = length(aName);
fwrite(aFid,aNameSize,'uint');
fwrite(aFid,aName,'char')';
