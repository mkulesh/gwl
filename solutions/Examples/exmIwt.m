
path(path, '../../mshell');

aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';

%---------------------------------------------------------------------------
figure(1) % Real signal - Shanon wavelet
aYmax = 0.5;
aDataType = 1;
gwlCreateAxis(1024,-3,3,'lin',aTimeName,'Time');
gwlCreateAxis(256,0.0001,6,'lin',aFreqName,'Frequency');

gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=shanon --wavpar=0.5 --time=0 --freq=2.5 --outtype=1 --outfile=' aSignalName])
[aTime,aSignal,aParSig] = gwlSignalRead(1,aSignalName,'func',['--format=ASCII --istime --mult=0.39 --nomess'],aSignalName,'Shanon wavelet');

gwlPlotFunction(aTime,aSignal,0.07,0.8,0.9,0.15,min(aTime),max(aTime),-aYmax,aYmax,' ',gwlGetNotation('SIG','T'),['(a) ' aParSig.aName ' (source signal)']);

gwlCwt(aDataType, aSignalName, aFreqName, 1, 'morlet', 5, aSpectrName);
[aTime,aFreq,aCwt,aParams] = gwlConvert('3,5','--filter=1',aSpectrName);

gwlPlotImage(aTime, aFreq, aCwt(:,:,1), 0.07,0.56,0.9,0.2, '', gwlGetNotation('FREQ'), ['(b) ' gwlGetNotation('SIG','WABS')]);
gwlPlotImage(aTime, aFreq, aCwt(:,:,2), 0.07,0.35,0.9,0.2, ' ', gwlGetNotation('FREQ'), ['(c) ' gwlGetNotation('SIG','WARG')]);

[aTimeInv,aSignalInv,aParInv] = gwlIwt(aDataType, aSpectrName, 'delta', 1, 2);

gwlPlotFunction(aTime,aSignal,0.07,0.15,0.9,0.15,min(aTime),max(aTime),-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('SIG','T'),['(d) ' aParInv.aName]);
    hold on; plot(aTimeInv,aSignalInv,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;

%---------------------------------------------------------------------------
figure(2); % Complex signal
aYmax = 10.0;
aDataType = 2;
gwlCreateAxis(1024,0,5.115,'lin',aTimeName,'Time');
gwlCreateAxis(256,0.0001,15,'lin --sign=full',aFreqName,'Frequency');

[aTime,aSignal,aParSig] = gwlSignalGen(aDataType,aTimeName,'harmrot','2.0,7.0,1.0,2.0,5.0,5.0,0.318',aSignalName,'Rotated harmonic function');
gwlPlotFunction(aTime,real(aSignal),0.07,0.75,0.9,0.15,min(aTime),max(aTime),-aYmax,aYmax,'',gwlGetNotation('MSIG','T','x'),['(a) ' aParSig.aName ' (source signal)']);
gwlPlotFunction(aTime,imag(aSignal),0.07,0.6,0.9,0.15,min(aTime),max(aTime),-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('MSIG','T','z'),'');

gwlCwt(aDataType, aSignalName, aFreqName, 1, 'morlet', 5, aSpectrName);
[aTimeInv,aSignalInv] = gwlIwt(aDataType, aSpectrName, 'delta');

gwlPlotFunction(aTime,real(aSignal),0.07,0.30,0.9,0.15,min(aTime),max(aTime),-aYmax,aYmax,'',gwlGetNotation('MSIG','T','x'),['(b) ' aParInv.aName]);
    hold on; plot(aTimeInv,real(aSignalInv),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;
gwlPlotFunction(aTime,imag(aSignal),0.07,0.15,0.9,0.15,min(aTime),max(aTime),-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('MSIG','T','z'),'');
    hold on; plot(aTimeInv,imag(aSignalInv),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;
      
%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete(aSignalName); delete(aSpectrName); 
clear all;

