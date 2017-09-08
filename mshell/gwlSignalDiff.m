function aDest = gwlSignalDiff(aSource, aSampleFreq)
% aDest = gwlSignalDiff(aSource, aSampleFreq)
% This is a simple procedure to differentiate a signal aSource having the sampling
% frequency aSampleFreq. Differentiation is performed along the first index of the
% source signal for each channel (second index) using a second order formula.
%
% Input parameters: 
%   aSource[m,n] is the source signal with m points and n channels
%   aSampleFreq = 1.0/(Time(2)-Time(1)) is the sampling frequency of the signal
%
% Output parameters:   
%   aDest[m,n] is the differentiated signal
% 
% Examples: 
%   To read and differentiate all columns of a file, the following code can be used:
%   [aTime,aSignal,aSigPar] = gwlSignalRead(1,'signal.asc','seis','--format=ASCII --istime');
%   aDiff = gwlSignalDiff(aSignal, aSigPar.aSample);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

[aCount,aChan]=size(aSource);
for n=1:aChan
    A1(1,n)=0;
    for k=2:aCount-1
        A1(k,n)=(aSampleFreq/2.0)*aSource(k+1,n)-(aSampleFreq/2.0)*aSource(k-1,n);
    end;
    A1(aCount,n)=A1(aCount-1,n);
    A1(1,n)=A1(2,n);
end;

aDest = A1;
