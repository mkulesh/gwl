function [aAxis,aSignal,aParams] = gwlReadSAC3D(aFiles)
% [aAxis,aSignal,aParams] = gwlReadSAC3D(aFiles)
% This procedure reads a 3-c signal file with format 'SAC'.
%
% Input parameters: 
%   aFiles is a structure with following fields:
%   aFiles.count=3; % We have three channels
%   aFiles.file1='chN.sac'; % Here we define the files
%   aFiles.file2='chE.sac';
%   aFiles.file3='chZ.sac';
%
% Output parameters:   
%   aAxis is an array contained time values
%   aSignal is a matrix contained signal values 
%   aParams contains the technical details about the aSignal variable
%
% Examples: To read a SAC signal file, use the following code
%   InpitFiles.count=3; % We have three channels
%   InpitFiles.file1='BHN.sac'; % Here we define the files as N-E-Z
%   InpitFiles.file2='BHE.sac';
%   InpitFiles.file3='BHZ.sac';
%   [aTime,DATA_RAW,aPar] = gwlReadSAC3D(aInpitFiles);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aSignal=[];

for  kfile=1:aFiles.count      
    eval(['filename = aFiles.file'  num2str(kfile) ]);
    [fid,message] = fopen(filename);
    [fheader,count] = fread(fid,70,'float');
    [nheader,count] = fread(fid,35,'long');
    [lheader,count] = fread(fid,5,'long');
    [kheader,count] = fread(fid,[8,24],'char');
    
    npts = nheader(10);
    b = fheader(6);
    if abs(b) > 0 
        'check if files are synchronized'
    end  
    e = fheader(7);
    deltaf(kfile) = fheader(1); % sample interval in sec
    if kfile==1
        delt=deltaf(1);
    else 
        if mean(deltaf)==deltaf(1)
            delt=deltaf(1); 
        else 
            'the files have different sample rate'
            kfile=100
        end 
    end 
    
    aParams.sr =1/delt; 
    date = nheader(1:2);
    
    aParams.yr(kfile)=date(1); 
    aParams.dy(kfile)=date(2);
    aParams.hour(kfile) = nheader(3);
    aParams.minu(kfile) = nheader(4);
    aParams.sec(kfile) = nheader(5) + (nheader(6)/1000);
    aParams.station = setstr(kheader(:,1)');
    aParams.stat = aParams.station(1:4);
    aParams.chan= aParams.station(5:8); 
    aParams.stla=fheader(32); 
    aParams.stlo=fheader(33); 
    aParams.evla= fheader(36); 
    aParams.evlo=fheader(37); 
    aParams.epicdist_deg=fheader(54); 
    aParams.otime=fheader(8); % event origin time (seconds relative to ref time);
    
    data = zeros(npts,1);
    [data,count] = fread(fid,npts,'float');
    data =data - mean(data);
    close_stat = fclose(fid);
    gw_a=data-mean(data); 
    
    [pD,qD]=size(aSignal);
    [pd,qd]=size(data); 
    if pD==0
        aSignal=[aSignal gw_a];
    else 
        if pD==pd
            aSignal=[aSignal gw_a];
        else 
            'data is adjusted to the minimum length to be the same length for each station'
            min_length=min(pd,pD); 
            aSignal=aSignal(1:min_length,:);
            gw_a=gw_a(1:min_length,:); 
            aSignal=[aSignal gw_a];
        end
    end
    
end;

[aParams.aPointCount,aParams.aChanCount] = size(aSignal);
aShift = 3600.*aParams.hour + 60.*aParams.minu + aParams.sec;
aParams.aStartSec = min(aShift);
aParams.aShift = floor(aParams.sr.*(aShift-aParams.aStartSec)) + 1;
aAxis = aParams.aStartSec + (0:(aParams.aPointCount-1))/aParams.sr;

clear gw_a;

