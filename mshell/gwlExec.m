function gwlExec(aProg,aPar)
% gwlExec(aProg,aPar)
% Run a GWL module from binary folder. The path to this folder is different for Windows and Linux versions.
%
% Input parameters: 
%   aProg is a string defines the name of GWL module
%   aPar is a string defines the parameters of this module
%
% Example:
%   gwlExec('gwlET3D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=acovar --filter=1 --tw=2']);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

str = computer;
if strcmp(str,'PCWIN') == 1
   aPath = '..\..\bin\';
   aExt = '.exe';
else
   aPath = '../../bin/';
   aExt = '';
end;    
aCommand = [aPath aProg aExt ' ' aPar ' --outtype=2'];
status = system(aCommand);

