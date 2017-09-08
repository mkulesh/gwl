function [aAxis,aSignal,aParams]=gwlIwt(aDataType, aSpectrName, aWavelet, aWavPar, aAmpl, aFile, aName)
% [aAxis,aSignal,aParams]=gwlIwt(aDataType, aSpectrName [, aWavelet, aWavPar, aAmpl, aFile, aName])
% Component-wise calculation of the inverse continuous wavelet transform corresponds to a multi-component 
% real or complex signal.
%
% Input parameters define source wavelet spectrum, algorithm, wavelet and its properties: 
%   aDataType defines the type of the output signal:
%     aDataType=1 corresponds to a real signal;
%     aDataType=2 - to a complex signal
%   aSpectrName is a file name with the source wavelet spectrum. This spectrum can be a multi-component spectrum.
%   aWavelet is a string variable that is the name of analyzing wavelet: 
%     'morlet', 'morletre', 'cauchy', 'shanon', 'haar' - correspond to a full time-frequency integral;
%     'delta' - corresponds to an fast reconstruction algorithm based on a frequency integral (by default).
%     It is also possible to use here additional parameter '--cutoff=C', where real parameter C defines 
%     the amplitude cut-off properties of the wavelet, by default --cutoff=0.01
%   aWavPar is a real corresponds to the time-frequency resolution of the wavelet
%   aFile is the string variable that defines the name of the output file for 'gwl' algorithm
%   aName is the name of the spectrum object that will be saved in the output file aFile
%   aAmpl is related to the normalization coefficient of the inverse transform.
%     We do not need to calculate this coefficient for all applications; if aAmpl is not given or
%     aAmpl=2, gwlIwt will calculate automatically and use this coefficient
%   aFile is the string variable that defines the name of the output file with reconstructed signal
%   aName is the name of the signal object that will be saved in the output file aFile
%
% Output parameters:   
%   aAxis is an array contained time axis of the wavelet spectrum 
%   aSignal is a complex or real array contains all channels of the reconstructed signal
%   aParams contains the technical details about the aSignal variable
%
% Examples:
%   1) Calculate a complex inverse wavelet transform using Morlet wavelet without normalization coefficient
%      [aTimeInv,aSignalInv,aParInv] = gwlIwt(2, 'spec.dat', 'morlet', 1, 1);
%   2) Calculate a real inverse wavelet transform using fast reconstruction formula with normalization coefficient
%      [aTimeInv,aSignalInv,aParInv] = gwlIwt(1, 'spec.dat', 'delta', 0, 2, 'signal.dat');
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aPar = '';

if (nargin == 7) 
    aPar = [' --name="' aName '"'];
end;

if (nargin < 6) 
    aFile = tempname;
end;
aPar = [aPar ' --outfile=' aFile];

if (nargin < 5)
    aPar = [aPar ' --ampl'];
else if (aAmpl == 2) 
    aPar = [aPar ' --ampl']; end;
end;

if (nargin > 3) 
     aPar = [aPar ' --wavpar=' num2str(aWavPar)];
end;

if (nargin > 2) 
     aPar = [aPar ' --wavelet=' aWavelet];
else
     aPar = [aPar ' --wavelet=delta'];
end;

if (aDataType==2)
    aPar = [aPar ' --iscmpl'];
end;

aPar = [aPar ' --infile=' aSpectrName];

gwlExec('gwlIwt',aPar);

if (nargout > 0) 
    fid = fopen(aFile,'r');
    if (nargout < 3) 
        [aAxis,aSignal]=gwlReadSignal(fid);
    else
        [aAxis,aSignal,aParams]=gwlReadSignal(fid);
    end;
    fclose(fid);
end;

if (nargin < 6) 
    delete(aFile);
end;
