function [aAxis,aParams]=gwlCreateAxis(aCount, aMin, aMax, aScale, aFile, aName)
% [aAxis,aParams]=gwlCreateAxis(aCount, aMin, aMax [, aScale, aFile, aName])
% Using this procedure, we can create a GWL axis object with help of bin/gwlCreateAxis program.
% The axis object is used in some modules of mshell library: in gwlCwt and 
% gwlDispModel as input frequency axis and in gwlSignalGen as input time axis.
%
% Input parameters define axis properties:
%   aCount - count of axis points
%   aMin - minimal axis value
%   aMax - maximal axis value
%   aScale is a string variable that can contain two parameters:
%      'lin' or 'log' define the linear or logarithmic scale of the axis;
%      '--sign' defines the symmetry of the axis and is important for gwlCwt procedure, 
%      where symmetry corresponds to spectrum type;
%      '--sign=full': axis will be symmetric (full spectrum);
%      '--sign=prog': axis will contain only positive values (prograde spectrum);
%      '--sign=reg': axis will contain only negative values (regrade spectrum)
%   aFile is the string variable that defines the name of the output file 
%   aName is the name of the axis object that will be saved in the output file aFile
%
% Output parameters:   
%   aAxis is an array contained axis values 
%   aParams contains the technical details about the aAxis variable
% 
% Examples: 
%   1) To produce a linear axis:
%      [aAxLin,aParLin] = gwlCreateAxis(1000,1,100,'lin','axislin.dat','Linear axis');
%      or
%      [aAxLin,aParLin] = gwlCreateAxis(1000,1,100);
%   2) To produce a logarithmic axis:   
%      [aAxLog,aParLog] = gwlCreateAxis(1000,1,100,'log','axislog.dat','Logarithmic axis');
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
    
if (nargin > 3) 
    aPar = [aPar ' --scale=' aScale];
end;

aPar = [aPar ' --count=' num2str(aCount) ' --min=' num2str(aMin) ' --max=' num2str(aMax)];

gwlExec('gwlCreateAxis',aPar);

if (nargout > 0) 
    fid = fopen(aFile,'r');
    if nargout == 1
        aAxis = gwlReadAxis(fid);
    else     
        [aAxis,aParams] = gwlReadAxis(fid);
    end;
    fclose(fid);
end;

if (nargin < 5) 
    delete(aFile);
end;
