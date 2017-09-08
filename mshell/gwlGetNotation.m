function aLabel = gwlGetNotation(aType, aPar1, aPar2, aPar3, aPar4)
% aLabel = gwlGetNotation(aType, aPar1, aPar2, aPar3, aPar4)
% This procedure returns a string variable aLabel that is assembled from
% primitive symbols, TEX-symbols and dimensions depends on aType parameter
%
% This procedure collects almost all notations used in examples from
% gwl/solutions directory. We use this procedure to make the labels of all 
% GWL plots uniform.
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aLabel = '';

if(strcmp(aType,'IND')==1)
    aLabel = 'Index';
end;

if(strcmp(aType,'TRN')==1)
    aLabel = 'Trace Number';
end;

if(strcmp(aType,'AVAL')==1)
    aLabel = 'Value';
end;

if(strcmp(aType,'TIME')==1)
    aStr0 = '';
    if (nargin == 1) aStr0 = ' (s)'; end;
    if (nargin > 1) aStr0 = [' (' aPar1 ')']; end;
    aLabel = ['Time' aStr0];
end;

if(strcmp(aType,'FREQ')==1)
    aLabel = 'Frequency (Hz)';
end;

if(strcmp(aType,'CFREQ')==1)
    aLabel = 'Frequency (rad/sec)';
end;

if(strcmp(aType,'DEPTH')==1)
    aLabel = 'kz/(2\pi)';
end;

if(strcmp(aType,'DISP')==1)
    aStr0 = '';
    aStr1 = ''; 
    if (nargin > 1) 
        if(strcmp(aPar1,'ATN')==1)  aStr0 = '\alpha';  end;
        if(strcmp(aPar1,'WN')==1)  aStr0 = 'k';  end;
        if(strcmp(aPar1,'IMWN')==1)  aStr0 = '\Im k';  end;
        if(strcmp(aPar1,'VEL')==1) aStr0 = 'Velocity';  end;
        if(strcmp(aPar1,'CP')==1)  aStr0 = 'C_p';  end;
        if(strcmp(aPar1,'CG')==1) aStr0 = 'C_g';  end;
        if(strcmp(aPar1,'CPN')==1)  aStr0 = 'C_p';  end;
        if(strcmp(aPar1,'CGN')==1) aStr0 = 'C_g';  end;
    end;
    if(nargin > 2) 
        if(strcmp(aPar2,'F')==1)  aStr1 = '(f)';  end;
        if(strcmp(aPar2,'CF')==1)  aStr1 = '(\omega)';  end;
        if(strcmp(aPar2,'CFN')==1)  aStr1 = '(\omega)/C_t';  end;
    end;
    if(nargin > 3) 
        aLabel = [aStr0 aStr1 aPar3];
    else
        aLabel = [aStr0 aStr1];
    end;
end;

if(strcmp(aType,'SIG')==1 | strcmp(aType,'CSIG')==1 | strcmp(aType,'MSIG')==1 | strcmp(aType,'TSIG')==1 | strcmp(aType,'VSIG')==1 | strcmp(aType,'ASIG')==1 | strcmp(aType,'GSIG')==1)
    aPref = ''; 
    aArg = ''; 
    aInd = '';
    aPow = '';
    aSig = '';
    if(nargin > 1) 
        if(strcmp(aPar1,'T')==1)  aArg = '(t)';  end;
        if(strcmp(aPar1,'Z')==1)  aArg = '(z)';  end;
        if(strcmp(aPar1,'RET')==1)  aPref = '\Re'; aArg = '(t)';  end;
        if(strcmp(aPar1,'IMT')==1)  aPref = '\Im'; aArg = '(t)';  end;
        if(strcmp(aPar1,'F')==1)  aPref = ''; aArg = '(f)';  end;
        if(strcmp(aPar1,'CF')==1)  aPref = ''; aArg = '(\omega)';  end;
        if(strcmp(aPar1,'REF')==1)  aPref = '\Re'; aArg = '(f)';  end;
        if(strcmp(aPar1,'IMF')==1)  aPref = '\Im'; aArg = '(f)';  end;
        if(strcmp(aPar1,'FM')==1)  aPref = '|'; aArg = '(f)|';  end;
        if(strcmp(aPar1,'WG')==1) aPref = 'W_g'; aArg = '(t,f)';  end;
        if(strcmp(aPar1,'WGP')==1) aPref = 'W^+_g'; aArg = '(t,f)';  end;
        if(strcmp(aPar1,'WGM')==1) aPref = 'W^-_g'; aArg = '(t,f)';  end;
        if(strcmp(aPar1,'WABS')==1) aPref = '|W_g'; aArg = '(t,f)|';  end;
        if(strcmp(aPar1,'WABSP')==1) aPref = '|W^+_g'; aArg = '(t,f)|';  end;
        if(strcmp(aPar1,'WABSM')==1) aPref = '|W^-_g'; aArg = '(t,f)|';  end;
        if(strcmp(aPar1,'WARG')==1) aPref = 'arg W_g'; aArg = '(t,f)';  end;
        if(strcmp(aPar1,'WARGP')==1) aPref = 'arg W^+_g'; aArg = '(t,f)';  end;
        if(strcmp(aPar1,'WARGM')==1) aPref = 'arg W^-_g'; aArg = '(t,f)';  end;
    end;
    if (nargin > 2 & strcmp(aPar2,'')==0) aInd = ['_{' num2str(aPar2) '}']; end;
    if (nargin > 3 & strcmp(aPar3,'')==0) aPow = ['^{' num2str(aPar3) '}']; end;
    if(strcmp(aType,'SIG')==1) aSig = 's'; end;
    if(strcmp(aType,'CSIG')==1) aSig = 'z'; end;
    if(strcmp(aType,'MSIG')==1) aSig = 'u'; end;
    if(strcmp(aType,'TSIG')==1) aSig = 'T'; end;
    if(strcmp(aType,'VSIG')==1) aSig = 'v'; end;
    if(strcmp(aType,'ASIG')==1) aSig = 'a'; end;
    if(strcmp(aType,'GSIG')==1) aSig = 'g'; end;
    if(nargin > 4)
        aLabel = [aPref aSig aPow aInd aArg aPar4]; 
    else
        aLabel = [aPref aSig aPow aInd aArg]; 
    end;
end;

if(strcmp(aType,'EPAR')==1)
    aStr0 = '';
    aStr1 = ''; 
    aInd = '';
    if (nargin > 1) 
        if(strcmp(aPar1,'RMIN')==1)   aStr0 = 'r';  end;
        if(strcmp(aPar1,'RMED')==1)   aStr0 = 'r_s';  end;
        if(strcmp(aPar1,'RMAX')==1)   aStr0 = 'R';  end;
        if(strcmp(aPar1,'RATIO')==1)  aStr0 = '\rho';  end;
        if(strcmp(aPar1,'RATIO1')==1) aStr0 = '\rho_1';  end;
        if(strcmp(aPar1,'RATIOS')==1) aStr0 = '\rho_s';  end;
        if(strcmp(aPar1,'TILT')==1)   aStr0 = '\theta';  end;
        if(strcmp(aPar1,'PDIFF')==1)  aStr0 = '\Delta\phi';  end;
        if(strcmp(aPar1,'DIP')==1)    aStr0 = '\beta';  end;
        if(strcmp(aPar1,'AZ')==1)     aStr0 = '\gamma';  end;
        if(strcmp(aPar1,'BAZ')==1)    aStr0 = '\gamma_0';  end;
        if(strcmp(aPar1,'ATIME')==1)  aStr0 = '\Deltat';  end;
        if(strcmp(aPar1,'INST')==1)   aStr0 = '\Omega';  end;
    end;
    if(nargin > 2) 
        if(strcmp(aPar2,'T')==1)  aStr1 = '(t)';  end;
        if(strcmp(aPar2,'F')==1)  aStr1 = '(f)';  end;
        if(strcmp(aPar2,'WG')==1) aStr1 = '(t,f)';  end;
    end;
    if (nargin > 3 & strcmp(aPar3,'')==0) aInd = ['_{' num2str(aPar3) '}']; end;
    if(nargin > 4)
        aLabel = [aStr0 aInd aStr1 aPar4];
    else
        aLabel = [aStr0 aInd aStr1];
    end;
end;

