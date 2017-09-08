function exmCwtMorlet()

%---------------------------------------------------------------------------
path(path, '../../mshell');
aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';
aDataType = 2; % 1=real, 2=complex
aFreq0 = 2; % frequency position of the wavelet (Hz)
gwlCreateAxis(1024,0,10,'lin',aTimeName,'Time');
gwlCreateAxis(256,0.001,10,'lin',aFreqName,'Frequency');

gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=morlet --wavpar=1 --time=5 --freq=' num2str(aFreq0) ' --outtype=1 --outfile=' aSignalName])
[aTime,aSignal,aPar] = gwlSignalRead(aDataType,aSignalName,'func',['--format=ASCII --istime --nomess'],aSignalName,'Morlet wavelet');
gwlPlotFunction(aTime,real(aSignal),0.07,0.69,0.9,0.25,min(aTime),max(aTime),-aFreq0,aFreq0,'',gwlGetNotation('CSIG','T'),'(a)',aPar.aName);
    hold on; plot(aTime,imag(aSignal),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;

gwlCwt(aDataType, aSignalName, aFreqName, 2, 'morlet', 0.5, aSpectrName);
[aTime,aFreq,aCwt,aParams] = gwlConvert('3,5','--filter=1',aSpectrName);

gwlPlotImage(aTime, aFreq, aCwt(:,:,1), 0.07,0.38,0.9,0.3, '', gwlGetNotation('FREQ'), ['(b) ' gwlGetNotation('CSIG','WABS')]);
gwlPlotImage(aTime, aFreq, aCwt(:,:,2), 0.07,0.07,0.9,0.3, gwlGetNotation('TIME'), gwlGetNotation('FREQ'), ['(c) ' gwlGetNotation('CSIG','WARG')]);

aParams

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete(aSignalName); delete(aSpectrName); 
clear all;

print -f1 -r600 -depsc exmCwtMorlet;
