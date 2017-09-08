function [aAxis,aSignal,aParams]=gwlSignalGen(aDataType, aAxisName, aType, aParams, aFile, aName)
% [aAxis,aSignal,aParams]=gwlSignalGen(aDataType, aAxisName, aType, aParams [, aFile, aName])
% This procedure generates a synthetic signal and writes it into a binary file compatible with GWL library.
%
% Input parameters define a type and parameters of generated synthetic signal:
%   aDataType defines the type of the output signal:
%     aDataType=1 corresponds to a real signal;
%     aDataType=2 - to a complex signal
%   aAxisName is the file name of the time axis used to generate the signal
%   aType corresponds to the type of synthetic signal with parameters aParams:
%     aType='zero' is a zero function;
%     aType='delta' is a delta function. aParams=t0,a0, where t0 is the position and a0 - the amplitude of the peak;
%     aType='harmon' is a real or complex harmonic function with amplitude change and phase difference:
%       F(t) = (p0+p1*e(t))*cos(2.0*PI*p2*t) + i(p3+p4*e(t))*sin(2.0*PI*p5*t+p6), e(t)=0..1, aParams=p0,p1,p2,p3,p4,p5,p6;
%     aType='harmphase' is a real or complex harmonic function with phase difference change:
%       F(t) = p0*sin(2.0*PI*p2*t+sin(p3*i)) + i*p1*sin(2.0*PI*p2*t), aParams=p0,p1,p2,p3;
%     aType='harmrot' is a real or complex rotated harmonic function, aParams=p0,p1,p2,p3,p4,p5,p6;
%     aType='rickdiss' is a real signal corresponds to a propagated Ricker wavelet, aParams=p0,p1,p2,p3,p4,p5
%   aFile is the string variable that defines the name of the output file 
%   aName is the name of the signal object that will be saved in the output file aFile
% 
% Output parameters:   
%   aAxis is the time axis related to the generated signal aSignal 
%   aParams contains the technical details about the aSignal variable
% 
% Examples: 
%   To generate a complex rotated harmonic function
%   [aTime,aSignal,aParSig] = gwlSignalGen(2,'time.dat','harmrot','2.0,7.0,1.0,2.0,5.0,5.0,0.318');
% 
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aPar = '';

if (nargin == 6) 
    aPar = [' --name="' aName '"'];
end;

if (nargin < 5) 
    aFile = tempname;
end;
aPar = [aPar ' --outfile=' aFile];

if (aDataType==2)
    aPar = [aPar ' --iscmpl'];
end;

aPar = [aPar ' --type=' aType ' --infile=' aAxisName ' --par=' aParams];

gwlExec('gwlSignalGen',aPar);

if (nargout > 0) 
    fid = fopen(aFile,'r');
    if (nargout < 3) 
        [aAxis,aSignal]=gwlReadSignal(fid);
    else
        [aAxis,aSignal,aParams]=gwlReadSignal(fid);
    end;
    fclose(fid);
end;

if (nargin < 5) 
    delete(aFile);
end;
