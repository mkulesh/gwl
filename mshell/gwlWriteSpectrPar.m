function gwlWriteSpectrPar(aFid,aParams)
% gwlWriteSpectrPar(aFid,aParams)
% This procedure writes a binary stream with format 'WP1.4'.
% In GWL do not exist files of 'WP1.4' format, such stream is only a part of
% other GWL binary files, for example, a part of spectrum container 'SP1.3'.
% Therefore, this procedure is only used inside of GWL mshell library.
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'wb') command.
%   aParams contains the technical details about the used wavelet, it is an
%   output of gwlReadSpectrPar procedure.
%
% This file is part of the GWL library. Copyright (C) 2006-2008 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'WP1.4';
aDataType = 8;
fwrite(aFid,aVer,'char')';
fwrite(aFid,aDataType,'uint');

fwrite(aFid,aParams.aTransType,'uint');
fwrite(aFid,aParams.aWavelet,'uint');
fwrite(aFid,aParams.aWaveletPar,'double');
fwrite(aFid,aParams.aIWavelet,'uint');
fwrite(aFid,aParams.aIWaveletPar,'double');
gwlWriteAxis(aFid,aParams.aFreq,aParams.aParFreq);
%fwrite(aFid,aParams.aTMin,'double');
%fwrite(aFid,aParams.aTMax,'double');
aNameSize = length(aParams.aWPName);
fwrite(aFid,aNameSize,'uint');
fwrite(aFid,aParams.aWPName,'char')';

