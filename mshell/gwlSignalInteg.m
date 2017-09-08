function aDest = gwlSignalInteg(aSource, aSampleFreq)
% aDest = gwlSignalInteg(aSource, aSampleFreq)
% This is a simple procedure to integrate a signal aSource having the sampling
% frequency aSampleFreq. Integration is performed along the first index of the
% source signal for each channel (second index) using a second order formula.
%
% Input parameters: 
%   aSource[m,n] is the source signal with m points and n channels
%   aSampleFreq = 1.0/(Time(2)-Time(1)) is the sampling frequency of the signal
%
% Output parameters:   
%   aDest[m,n] is the integrated signal
% 
% Examples: 
%   To read and integrate all columns of a file, the following code can be used:
%   [aTime,aSignal,aSigPar] = gwlSignalRead(1,'signal.asc','seis','--format=ASCII --istime');
%   aInteg = gwlSignalInteg(aSignal, aSigPar.aSample);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

[aCount,aChan]=size(aSource);
for n=1:aChan
    A1(1,n)=0;
    for k=2:aCount-1
        A1(k,n)=(A1(k-1,n)+(aSource(k+1,n) + 4.0*aSource(k,n) + aSource(k-1,n))/(6.0*aSampleFreq));
    end;
    A1(aCount,n)=0;
end;

aDest = A1;
