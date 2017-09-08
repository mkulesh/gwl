function [aFreq,aFour,aParams] = gwlFourTrans(aAlg, aSource, aShift, aSampleFreq, aFile, aName)
% [aFreq,aFour,aParams] = gwlFourTrans(aAlg, aSource, aShift [, aSampleFreq, aFile, aName])
% Calculation of the Fourier transform of a signal using both standard MATLAB fft()
% subroutine and gwlCft module.
%
% Input parameters define source signal, algorithm and its properties: 
%   aAlg is the string defined the algorithm of FT calculation: 
%     aAlg = 'mat' - the MATLAB fft() procedure is used;
%     aAlg = 'gwl' - the gwlCft module is used
%   aSource is the string defined the input signal. For MATLAB algorithm it is an array, for gwlCft it is a file name
%   aShift defines shift scheme:
%     aShift=1 - all frequencies are positive without shift;
%     aShift=2 - modify the spectrum to regrade/prograde form
%   aSampleFreq - sampling frequency of the source signal in the case of aAlg = 'mat'; for 'gwl' algorithm can be zero.
%   aFile is the string variable that defines the name of the output file for 'gwl' algorithm
%   aName is the name of the spectrum object that will be saved in the output file aFile
%
% Output parameters:   
%   aFreq is an array contained frequency values 
%   aFour is a complex matrix contains the Fourier coefficients for all channels of the signal aSource
%   aParams contains the technical details about the aFour variable
% 
% Examples: fft of a complex signal from the file 'signal.asc' in regrade/prograde form
%   aShift = 2;
%   aSignalName = 'signal.asc';
%   1) MATLAB fft:
%      [aTime,aSignal,aParSig] = gwlSignalRead(2,aSignalName,'func','--format=ASCII --istime'); 
%      [aFreq2,aFour2] = gwlFourTrans('mat',aSignal,aShift,aParSig.aSample);
%   2) gwlCft:
%      [aFreq1,aFour1] = gwlFourTrans('gwl',aSignalName,aShift);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

if aAlg == 'mat'
    aCount=length(aSource);
    if aShift == 1
        for k=1:aCount
            aFreq(k)=(aSampleFreq*(k-1)/(aCount-1));
        end;
        aFour = fft(aSource);
    else
        for k=1:aCount
            aFreq(k)=(aSampleFreq*(k-1)/(aCount-1) - aSampleFreq/2);
        end;
        aFour = fftshift(fft(aSource));
    end;
else
    aPar = ['--infile=' aSource];
    if (nargin == 6) 
        aPar = [' --name="' aName '"'];
    end;
    if (nargin < 5) 
        aFile = tempname;
    end;
    aPar = [aPar ' --outfile=' aFile];
    if aShift ~= 1
        aPar = [aPar ' --shift'];
    end;
    
    gwlExec('gwlCft',aPar);
    
    if (nargout > 0) 
        fid = fopen(aFile,'r');
        if (nargout < 3) 
            [aFreq,aFour]=gwlReadSignal(fid);
        else
            [aFreq,aFour,aParams]=gwlReadSignal(fid);
        end;
        fclose(fid);
    end;
    
    if (nargin < 5) 
        delete(aFile);
    end;
end;

