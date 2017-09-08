function [aFreq,aModel,aParams]=gwlDispModel(aFreq, aWn, aWnPar, aAtn, aAtnPar, aFile, aName)
% [aFreq,aModel,aParams]=gwlDispModel(aFreq, aWn, aWnPar [, aAtn, aAtnPar, aFile, aName])
% This function allows to create a dispersion model that is used in the module bin/gwlDiffeoDisp.
%
% Input parameters define the frequency-dependent wavenumber and attenuation functions: 
%   aFreqName is a file name with the frequency axis that will be use to produce wavenumber and attenuation.
%     This axis can be created using gwlCreateAxis procedure. See help of gwlCreateAxis for details.
%   aWn and aWnPar are string and real array variables that define the type of wavenumber approximation
%     and its parameters:
%     aWn = 'vel' is a Gaussian phase velocity approximation;
%     aWn = 'gauss' is a exponential approximation with 3 parameters;
%     aWn = 'polin' is a polynomial approximation with N parameters;
%     aWn = 'bspline' is a B-spline approximation with N parameters;
%     aWn = 'colecole' is a Cole-Cole model with 4 parameters
%   aAtn and aAtnPar are string and real array variables that define the type of wavenumber approximation
%     and its parameters:
%     aWn = 'gauss' is a exponential approximation with 3 parameters;
%     aAtn = 'polin' is a polynomial approximation with N parameters;
%     aAtn = 'bspline' is a B-spline approximation with N parameters;
%     aAtn = 'colecole' is a Cole-Cole model with 4 parameters.
%     Note that, in the case of independent wavenumber and attenuation such as for exponential, 
%     polynomial and B-spline approximations, we can combine different types of aWn and aAtn
%   aFile is the string variable that defines the name of the output file for 'gwl' algorithm
%   aName is the name of the spectrum object that will be saved in the output file aFile
%
% Output parameters:   
%   aFreq is an array contained frequency values 
%   aModel is a matrix corresponds to six frequency-dependent dispersion curves:
%      aModel(:,1) is the wave number;
%      aModel(:,2) is the wave number derivative;
%      aModel(:,3) is the phase velocity;
%      aModel(:,4) is the group velocity;
%      aModel(:,5) is the attenuation;
%      aModel(:,6) is the attenuation derivative
%   aParams contains the technical details about the aModel variable
%
% Examples: 
%   1) To produce a dispersion model with Gaussian phase velocity approximation and constant attenuation
%      [aFreq, aModel] = gwlDispModel('freq.dat', 'vel', '1300,400,30', 'polin', '0', 'model.dat');
%   2) To produce a dispersion model with Gaussian phase velocity approximation and polynomial attenuation
%      [aFreq, aMode1] = gwlDispModel('freq.dat', 'vel', '1256,-62,50', 'polin', '6.204e-04,3.698e-05,-3.200e-06,6.000e-08,-3.336e-10', 'mode1.dat');
%   3) To produce a dispersion Cole-Cole model
%      [aFreq, aModel] = gwlDispModel('freq.dat', 'colecole', '7.87E+06,0.4,4.73E-04,1.717E-04', 'colecole', '0', 'model.dat');
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de



aPar = '';

if (nargin == 7) 
    aPar = [' --name="' aName '"'];
end;

if (nargin < 6) 
    aFile = tempname;
end;
aPar = [aPar ' --outfile=' aFile];
    
if (nargin > 3) 
    aPar = [aPar ' --atn=' aAtn ' --atnpar=' aAtnPar];
end;

aPar = [aPar ' --infile=' aFreq ' --wn=' aWn ' --wnpar=' aWnPar];

gwlExec('gwlDispModel',aPar);

if (nargout > 0) 
    fid = fopen(aFile,'r');
    if nargout == 2
        [aFreq, aModel] = gwlReadDispModel(fid);
    else     
        [aFreq, aModel, aParams] = gwlReadDispModel(fid);
    end;
    fclose(fid);
end;

if (nargin < 6) 
    delete(aFile);
end;
