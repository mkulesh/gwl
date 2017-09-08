function evalOneStation(aFileName,aChanNot,aChann,aTimeInt,aFreqName,aNoConvert)
% evalOneStation(): Polarization calculation subroutine used in eventSeptember20.m, filterHER.m, twoStations.m and exmCwtHER.m.

aNotation = strcat(aFileName,aChanNot);
aSignalName = strcat(aNotation,'sig.dat');
aElliparName = strcat(aNotation,'par.dat');
aSpectrName = strcat(aNotation,'cwt.dat');;
[aTime,aSigMSR] = gwlSignalRead(2,strcat(aFileName,'.asc'),'func',['--format=ASCII --chan=' aChann ' --tmin=' num2str(aTimeInt(1)) ' --tmax=' num2str(aTimeInt(2)) ' --istime'],aSignalName,aNotation);
gwlCwt(2, aSignalName, aFreqName, 1, 'cauchy', 16, aSpectrName,aNotation);
gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=complex --filter=15 --name=' aNotation]);
if (nargin == 5) 
    gwlConvert('13,4,14','--degree',aElliparName,aElliparName,aNotation);
end;

