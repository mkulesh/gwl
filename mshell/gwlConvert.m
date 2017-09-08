function [aTime,aFreq,aCwt,aParams]=gwlConvert(aComp, aConvpar, aSpectrName, aFile, aName)
% [aTime,aFreq,aCwt,aParams]=gwlConvert(aComp, aConvpar, aSpectrName [, aFile, aName]):
% This procedure converts a complex or elliptical spectrum to a real representation.
%
% Input parameters: 
%   aSpectrName is the string variable that defines the name of the source spectrum
%   aComp='n1,n2,...nn' is the string variable that contains comma-separates numbers. Each of this number 
%      defines a transformation that is related to the source format:
%      1) aSpectrName is the name of the COMPLEX object (for example, wavelet spectrum) 
%         n = 1:  real part
%         n = 2:  imaginary part
%         n = 3:  modulus
%         n = 4:  argument ('--degree' in aConvpar converts it from radians to degrees)
%         n = 5:  argument with cut-off ('--filter=p' in aConvpar defines the cut-off range in percent to the 
%                 max. modulus, '--degree' in aConvpar converts it from radians to degrees)
%         n = 6:  normalized modulus
%         n = 7:  modulated argument
%      2) aSpectrName is the name of the 2-D POLARIZATION object
%         n = 1:  instantaneous major axis
%         n = 2:  instantaneous minor axis
%         n = 4:  ellipticity ratio
%         n = 7:  instantaneous phase x
%         n = 8:  instantaneous phase y
%         n = 13: phase difference ('--degree' in aConvpar converts it from radians to degrees)
%         n = 14: tilt angle ('--degree' in aConvpar converts it from radians to degrees)
%      3) aSpectrName is the name of the 3-D POLARIZATION object
%         n = 1:  instantaneous major axis
%         n = 2:  instantaneous minor axis
%         n = 3:  instantaneous second minor axis
%         n = 4:  ellipticity ratio
%         n = 5:  second ellipticity ratio: second minor/minor
%         n = 6:  third ellipticity ratio: second minor/major
%         n = 7:  instantaneous phase x
%         n = 8:  instantaneous phase y
%         n = 9:  instantaneous phase z
%         n = 10: planarity angle x-component
%         n = 11: planarity angle y-component
%         n = 12: planarity angle z-component
%         n = 15: directional cosine: planarity-X angle
%         n = 16: directional cosine: planarity-Y angle
%         n = 17: directional cosine: planarity-Z angle
%         n = 18: sign of ratio
%         n = 19: dip angle of rmax ('--degree' in aConvpar converts it from radians to degrees)
%         n = 20: azimuth of rmax ('--degree' in aConvpar converts it from radians to degrees,
%                 '--modpi' in aConvpar applies modulo pi operator to azimuth) 
%  aConvpar is a string that contains additional parameters of the modulus bin/gwlConvert:
%     '--degree' converts the output from radians to degrees (valid only for complex argument, phase 
%        difference, tilt angle, dip angle and azimuth);
%     '--filter=p' defines the cut-off range of complex argument in percent to the maximum modulus;
%     '--modpi' applies modulo pi operator to azimuth;
%     '--chan=c' get only one channel c of the source spectrum;
%     '--voices=v' defines that only v voices will be written in the output spectrum;
%     '--points=p' defines that only p points will be written in the output spectrum
%   aFile is the string variable that defines the name of the output file 
%   aName is the name of the spectrum object that will be saved in the output file aFile
%
% Output parameters:   
%   aTime and aFreq are the time and frequency axes related to the spectrum aCwt 
%   aParams contains the technical details about the aCwt variable
%
% Examples: 
% 1) Get modulus and argument in degrees from a complex wavelet spectrum
%    [aTime,aFreq,aCwt] = gwlConvert('3,4','--degree','cwt.dat');
% 2) Get modulus and argument with cut-off (2% of maximum modulus) from a complex wavelet spectrum
%    [aTime,aFreq,aCwt] = gwlConvert('3,5','--filter=2','cwt.dat');
% 3) Get ellipticity ratio, phase difference and tilt angle (in degrees) from a 2D polarization spectrum
%    [aTime,aFreq,aElli] = gwlConvert('4,13,14','--degree','2Dellipar.dat');
% 4) Get ellipticity ratio, third ellipticity ratio, dip angle and azimuth (in degrees) from a 3D 
%    polarization spectrum. Azimuth is converted to modulo pi. 
%    [aTime,aFreq,aElli] = gwlConvert('4,6,19,20','--degree --modpi','3Dellipar.dat');
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aPar = '';

if (nargin == 5) 
    aPar = [' --name="' aName '"'];
end;

if (nargin < 4) 
    aFile = tempname;
end;
aPar = [aPar ' --outfile=' aFile];

if (isnumeric(aConvpar) == 1) % for compatibility with ver. 1.3-1.4
    aPar = [aPar ' --infile=' aSpectrName ' --comp=' aComp ' --filter=' num2str(aConvpar)];
else
    aPar = [aPar ' --infile=' aSpectrName ' --comp=' aComp ' ' aConvpar];
end;

gwlExec('gwlConvert',aPar);

if (nargout > 0) 
    fid = fopen(aFile,'r');
    if (nargout < 4) 
        [aTime,aFreq,aCwt]=gwlReadSpectrum(fid);
    else
        [aTime,aFreq,aCwt,aParams]=gwlReadSpectrum(fid);
    end;
    fclose(fid);
end;

if (nargin < 4) 
    delete(aFile);
end;