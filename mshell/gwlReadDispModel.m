function [aAxis,aSignal,aParams]=gwlReadDispModel(aFid)
% [aAxis,aSignal,aParams]=gwlReadDispModel(aFid)
% This procedure reads a stream with format 'DM1.3' from a binary data file of GWL .
% The file with this format is produces in the module 'bin/gwlDispModel'
%
% Input parameters: 
%   aFid - a pointer on the data file, initialized before using aFid=fopen(aFile,'r') command.
%
% Output parameters:   
%   aAxis is an array contained frequencies values for dispersion curves
%   aSignal is a matrix corresponds to six frequency-dependent dispersion curves:
%      aSignal(:,1) is the wave number;
%      aSignal(:,2) is the wave number derivative;
%      aSignal(:,3) is the phase velocity;
%      aSignal(:,4) is the group velocity;
%      aSignal(:,5) is the attenuation;
%      aSignal(:,6) is the attenuation derivative
%   aParams contains the technical details about the aSignal variable
%
% Examples: 
%   To produce and read a dispersion model, use the following code
%   gwlExec('gwlDispModel', '--infile=freq.dat --outfile=disp.dat --analyt --wn=vel --wnpar=1100,400,20 --name="1th dispersion model"');
%   fid = fopen('disp.dat','r'); 
%   [aFreq,aModel]=gwlReadDispModel(fid); 
%   fclose(fid);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aVer = 'DM1.3';
aPar.aObjVer = fread(aFid,5,'*char')';
aPar.aDataType = fread(aFid,1,'uint');
if strcmp(aPar.aObjVer,aVer)==0
    fclose(aFid);
    error(['version or format of input file is incorrect in procedure: gwlReadDispModel(). Procedure needs ver. ' aVer]);
end;

aPar.aAnalytical = fread(aFid,1,'uint');
if aPar.aAnalytical==1
    aPar.aWn = fread(aFid,1,'uint');
    aPar.aAtn = fread(aFid,1,'uint');
    aPar.aWnPar = gwlReadVector(aFid);
    aPar.aWnAtn = gwlReadVector(aFid);
end;
aPar.aFreqMin = fread(aFid,1,'*double');
aPar.aFreqMax = fread(aFid,1,'*double');
[aAxis,aSignal]=gwlReadSignal(aFid);

if nargout == 3
    aParams = aPar;
end;
