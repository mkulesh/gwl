function [aTime,aFreq,aCwt,aParams]=gwlCwt(aDataType, aSignalName, aFreqName, aAlg, aWavelet, aWavPar, aFile, aName)
% [aTime,aFreq,aCwt,aParams]=gwlCwt(aDataType, aSignalName, aFreqName, aAlg [, aWavelet, aWavPar, aFile, aName])
% Calculation of the continuous wavelet transform of a multi-component real or complex signal component-wise.
%
% Input parameters define source signal, algorithm, analyzing wavelet and its properties: 
%   aDataType defines the type of the input signal:
%     aDataType=1 corresponds to a real signal;
%     aDataType=2 - to a complex signal
%   aSignalName is a file name with the source signal; this signal can be a multi-component seismogram
%   aFreqName is a file name with the frequency axis, this axis can be created using gwlCreateAxis procedure. 
%     The symmetry property of the axis defines that the wavelet spectrum will be prograde, 
%     regrade or full. See help of gwlCreateAxis for details
%   aAlg is an integer that defines the algorithm of CWT calculation:
%     aAlg=0 corresponds the slow, but precise integration in the time domain; 
%     aAlg=1 corresponds the FFT-based without automatic adjustment of signal and wavelet length (with edge reflection);
%     aAlg=2 corresponds the FFT-based with automatic adjustment of signal and wavelet length (without edge reflection)
%   aWavelet is a string variable that is the name of analyzing wavelet: 
%     'morlet', 'morletre', 'cauchy', 'shanon', 'haar'; by default aWavelet='morlet'.
%     It is also possible to use here additional parameter '--cutoff=C', where real parameter C defines 
%     the amplitude cut-off properties of the wavelet. By default --cutoff=0.01
%   aWavPar is a real corresponds to the time-frequency resolution of the wavelet
%   aFile is the string variable that defines the name of the output file with wavelet spectrum
%   aName is the name of the spectrum object that will be saved in the output file aFile
%
% Output parameters:   
%   aTime is an array contained time axis of the wavelet spectrum 
%   aFreq is an array contained frequency values 
%   aCwt is a complex array of matrices contain the wavelet coefficients for all channels of the source signal
%   aParams contains the technical details about the aCwt variable
%
% Examples:
%   1) Calculate the wavelet spectrum of a real signal from the file 'signal.dat' using Morlet wavelet with 
%      the parameter 0.5, linear frequency axis and FFT-based algorithm into the file 'spectr.dat':
%      gwlCreateAxis(256,0.001,10,'lin','freq.dat');
%      gwlCwt(1, 'signal.dat', 'freq.dat', 2, 'morlet', 0.5, 'spectr.dat');
%   2) Calculate the full wavelet spectrum of a complex signal from the file 'signal.dat' using Cauchy wavelet with 
%      the parameter 20 without cut-off, logariphmic frequency axis and FFT-based algorithm into the file 'spectr.dat':
%      gwlCreateAxis(256,0.001,10,'log --sign=full','freq.dat');
%      gwlCwt(2, 'signal.dat', 'freq.dat', 2, 'cauchy --cutoff=0', 20, 'spectr.dat');
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aPar = '';

if (nargin == 8) 
    aPar = [' --name="' aName '"'];
end;

if (nargin < 7) 
    aFile = tempname;
end;
aPar = [aPar ' --outfile=' aFile];

if (nargin > 5) 
     aPar = [aPar ' --wavpar=' num2str(aWavPar)];
end;

if (nargin > 4) 
     aPar = [aPar ' --wavelet=' aWavelet];
end;

if (aDataType==2)
    aPar = [aPar ' --iscmpl'];
end;

aPar = [aPar ' --infile=' aSignalName ' --freq=' aFreqName ' --wttype=' num2str(aAlg)];

gwlExec('gwlCwt',aPar);

if (nargout > 0) 
    fid = fopen(aFile,'r');
    if (nargout < 4) 
        [aTime,aFreq,aCwt]=gwlReadSpectrum(fid);
    else
        [aTime,aFreq,aCwt,aParams]=gwlReadSpectrum(fid);
    end;
    fclose(fid);
end;

if (nargin < 7) 
    delete(aFile);
end;
