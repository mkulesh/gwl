function [aAxis,aSignal,aParams]=gwlSignalRead(aDataType, aSourceName, aType, aParams, aFile, aName)
% [aAxis,aSignal,aParams]=gwlSignalRead(aDataType, aSourceName, aType, aParams [, aFile, aName]): 
% This procedure reads an ASCII file with the source signal and converts it into 
% an internal binary format of GWL library. Various formats of the input ASCII 
% file may be read with gwlSignalRead.
% 
% Input parameters: 
%   aDataType defines the type of output signal aSignal: aDataType=1 corresponds to 
%     the real-valued signal and aDataType=2 defines the complex signal
%   aSourceName is the string variable that defines the name of the source ASCII file
%   aType corresponds to the format of the output signal: aType=func defines aSignal 
%     to be a function (real or complex), aType=seis defines aSignal to be a 
%     seismogram and aType=seis2D corresponds to a 2-component seismogram.
%   aParams is a string that contains the additional command line parameters of the 
%     module bin/gwlSignalRead. Important are the following parameters: 
%     --format defines the input data format. We recommend using only default value 
%       --format=ASCII 
%     --istime indicates that the time column is presented in the source file as FIRST 
%       column
%     --smplfreq=S has to be used if the time column absents in the source file and 
%       --smplfreq defines in this case the sampling frequency of the signal; 
%       for example --smplfreq=100
%     --tmin=T1 and --tmax=T2 define the time interval. The signal only from the 
%       interval [T1...T2] will be read from the source file
%     --to2p indicates that the output signal will be zero-padded to the nearest 
%       2^n length. It is important for FFT-based processing
%     --chan=c1,c2,...,cn defines the sequence of the columns in the source file, 
%       which will be read into the output signal aSignal 
%     --mult=m is a coefficient to be multiplied with all channels of the output signal 
%     --resample=R allows us to resample the signal by R times
%     --rot=c1,c2,alpha performs a rotation of channels c1 and c2 by the angle alpha.
%     --nomess indicates that no messages will be showed 
%   aFile is the string variable that defines the name of the output file 
%   aName is the name of the signal object that will be saved in the output file aFile
% 
% Output parameters:   
%   aAxis is the time axis related to the output signal aSignal 
%   aParams contains the technical details about the aSignal variable
% 
% Examples: 
% 1) Read a two columns file as a function. First column is the time. Output 
%    signal will be multiplied with 0.1. 
%    [aTime,aSignal] = gwlSignalRead(1,'signal.asc','func','--format=ASCII --istime --mult=0.1 --nomess'); 
% 2) Read all columns of a file as seismogram. First column is the time. 
%    [aTime,aSignal] = gwlSignalRead(1,'signal.asc','seis','--format=ASCII --istime','signal.dat','Experimental signal'); 
% 3) Read all columns of a file as seismogram. No time is presented in the source file. 
%    [aTime,aSeis] = gwlSignalRead(1,'signal.asc','seis','--format=ASCII --smplfreq=4000'); 
% 4) Read four columns from a file as seismogram. First column is the time. Output 
%    signal has 2^n length. Only the data between 150 and 250 sec will be read from 
%    the source file. The output signal at t>250 sec is zero-padded. 
%    [aTime,aSignal] = gwlSignalRead(1,'signal.asc','seis','--chan=2,1,0 --tmin=150 --tmax=250 --to2p --istime');
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

aPar = [aPar ' --type=' aType ' --infile=' aSourceName ' ' aParams];

gwlExec('gwlSignalRead',aPar);

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
