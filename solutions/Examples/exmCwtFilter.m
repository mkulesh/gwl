function exmCwtFilter()

%---------------------------------------------------------------------------
path(path, '../../mshell');
aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';

figure(1); %---------------------------------------------------------------------------
% source signal - real (aDataType=1) Morlet wavelet with the frequency aFreq0
aDataType = 1; % 1=real, 2=complex
aFreq0 = 2; % frequency position of the wavelet (Hz)
gwlCreateAxis(1024,0,10,'lin',aTimeName,'time');
gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=morlet --wavpar=1 --time=5 --freq=' num2str(aFreq0) ' --outtype=1 --outfile=' aSignalName])
[aTime,aSignal] = gwlSignalRead(aDataType,aSignalName,'func',['--format=ASCII --istime'],aSignalName,'Morlet wavelet');
gwlPlotFunction(aTime,aSignal,0.07,0.69,0.9,0.25,min(aTime),max(aTime),-aFreq0,aFreq0,'','','Soure signal');

% wavelet spectrum of the source fignal
gwlCreateAxis(256,0.001,10,'lin',aFreqName,'frequency');
gwlCwt(aDataType, aSignalName, aFreqName, 1, 'morlet', 0.5, aSpectrName);
fid = fopen(aSpectrName,'rb'); 
    [aTime,aFreq,aCwt,aCwtParams] = gwlReadSpectrum(fid); 
    fclose(fid);
gwlPlotImage(aTime, aFreq, abs(aCwt), 0.07,0.38,0.9,0.3, '', 'Frequency', 'Argument');
gwlPlotImage(aTime, aFreq, angle(aCwt), 0.07,0.07,0.9,0.3, 'Time', 'Frequency', 'Phase');
aCwtParams

figure(2); %---------------------------------------------------------------------------
% now we calculate a filter for the wavelet spectrum and write it as a binary file
aAmax = max(max(abs(aCwt)));
for v=1:aCwtParams.aVoices
    for p=1:aCwtParams.aPoints
        if abs(aCwt(v,p))>0.9*aAmax
            aCwt(v,p) = 0;
        end;
    end;
end;
gwlPlotImage(aTime, aFreq, abs(aCwt), 0.07,0.64,0.9,0.3, '', 'Frequency', 'Argument');
gwlPlotImage(aTime, aFreq, angle(aCwt), 0.07,0.33,0.9,0.3, '', 'Frequency', 'Phase');
fid = fopen(aSpectrName,'wb'); 
    gwlWriteSpectrum(fid,aTime,aFreq,aCwt,aCwtParams); 
    fclose(fid);
    
% inverse wavelet transform    
[aTimeInv,aSignalInv] = gwlIwt(aDataType, aSpectrName, 'delta', 1, 2);
gwlPlotFunction(aTimeInv,aSignalInv,0.07,0.07,0.9,0.25,min(aTimeInv),max(aTimeInv),-aFreq0,aFreq0,'Time','','Filtered signal');
    hold on;   
    plot(aTime,aSignal,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    
    hold off;

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete(aSignalName); delete(aSpectrName); 
clear all;
