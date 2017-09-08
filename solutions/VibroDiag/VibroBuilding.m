function VibroBuilding()
% VibroBuilding(): Plotting of experimental data: spectral properties of 
% velocities and accelerations by a vibrational measurement at a building

%---------------------------------------------------------------------------
path(path, '../../mshell');
[aTime,aSignal,aPar] = gwlSignalRead(1,'VibroBuilding.asc','seis','--format=ASCII --smplfreq=100 --tmin=0 --tmax=327.67 --to2p');
aSFreq = aPar.aSample;

%---------------------------------------------------------------------------
figure(2);
aX1 = aSignal(:,1)/826.0;
aX2 = aSignal(:,2)/795.3;
aX3 = aSignal(:,3)/796.0;
aTmin = 20.0;
aTmax = 110.0;
aYmax = 0.0015;
gwlPlotFunction(aTime,aX1,0.07,0.67,0.9,0.25,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('VSIG','T','s','',', (m/s)'),'(a) Geophone: 1, Channel: S');
gwlPlotFunction(aTime,aX2,0.07,0.37,0.9,0.25,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('VSIG','T','w','',', (m/s)'),'(b) Geophone: 1, Channel: W');
gwlPlotFunction(aTime,aX3,0.07,0.07,0.9,0.25,aTmin,aTmax,-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('VSIG','T','z','',', (m/s)'),'(c) Geophone: 1, Channel: Z');

%---------------------------------------------------------------------------
figure(3);
aX1 = aSignal(:,4)/826.0;
aX2 = aSignal(:,5)/795.3;
aX3 = aSignal(:,6)/796.0;
aTmin = 0.0;
aTmax = 320.0;
aYmax = 0.01;
[aFreq,aF1] = gwlFourTrans('mat',aX1,2,aSFreq);
[aFreq,aF2] = gwlFourTrans('mat',aX2,2,aSFreq);
[aFreq,aF3] = gwlFourTrans('mat',aX3,2,aSFreq);
gwlPlotFunction(aTime,aX1,0.07,0.67,0.58,0.24,aTmin,aTmax,-aYmax,aYmax,' ',gwlGetNotation('VSIG','T','s','',', (m/s)'),'(a) Geophone: 4, Channel: S');
    gwlPlotFunction(aFreq,abs(aF1),0.7,0.67,0.25,0.24,0.0,50,0,max(abs(aF1)),' ','F_v');
gwlPlotFunction(aTime,aX2,0.07,0.37,0.58,0.24,aTmin,aTmax,-aYmax,aYmax,' ',gwlGetNotation('VSIG','T','w','',', (m/s)'),'(b) Geophone: 4, Channel: W');
    gwlPlotFunction(aFreq,abs(aF2),0.7,0.37,0.25,0.24,0.0,50,0,max(abs(aF2)),' ','F_v');
gwlPlotFunction(aTime,aX3,0.07,0.07,0.58,0.24,aTmin,aTmax,-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('VSIG','T','z','',', (m/s)'),'(c) Geophone: 4, Channel: Z');
    gwlPlotFunction(aFreq,abs(aF3),0.7,0.07,0.25,0.24,0.0,50,0,max(abs(aF3)),gwlGetNotation('FREQ'),'F_v');

%---------------------------------------------------------------------------
figure(4);
aYmax = 0.4;
aA1 = gwlSignalDiff(aX1,aSFreq);   [aFreq,aF1] = gwlFourTrans('mat',aA1,2,aSFreq);
aA2 = gwlSignalDiff(aX2,aSFreq);   [aFreq,aF2] = gwlFourTrans('mat',aA2,2,aSFreq);
aA3 = gwlSignalDiff(aX3,aSFreq);   [aFreq,aF3] = gwlFourTrans('mat',aA3,2,aSFreq);
gwlPlotFunction(aTime,aA1,0.07,0.67,0.58,0.24,aTmin,aTmax,-aYmax,aYmax,' ',gwlGetNotation('ASIG','T','s','',', (m/s^2)'),'(d) Geophone: 4, Channel: S');
    gwlPlotFunction(aFreq,abs(aF1),0.7,0.67,0.25,0.24,0.0,50,0,max(abs(aF1)),' ','F_a');
gwlPlotFunction(aTime,aA2,0.07,0.37,0.58,0.24,aTmin,aTmax,-aYmax,aYmax,' ',gwlGetNotation('ASIG','T','w','',', (m/s^2)'),'(e) Geophone: 4, Channel: W');
    gwlPlotFunction(aFreq,abs(aF2),0.7,0.37,0.25,0.24,0.0,50,0,max(abs(aF2)),' ','F_a');
gwlPlotFunction(aTime,aA3,0.07,0.07,0.58,0.24,aTmin,aTmax,-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('ASIG','T','z','',', (m/s^2)'),'(f) Geophone: 4, Channel: Z');
    gwlPlotFunction(aFreq,abs(aF3),0.7,0.07,0.25,0.24,0.0,50,0,max(abs(aF3)),gwlGetNotation('FREQ'),'F_a');

%---------------------------------------------------------------------------
figure(5);
aSource.Name = 'VibroBuilding.velD1z';      % source file name
aSource.Tmin = 500;                         % start time (in sec)
aSource.Tmax = 540;                         % end time (in sec)
aSource.Bin = 'signal.dat';                 % name of the temporary file with analysing signal
aWav.aFName = 'freq.dat';                   % File with frequency axis
aWav.aFCount = 256;                         % How many frequency points
aWav.aFMin = 0.0001;                        % Minimal frequency in Hz
aWav.aFMax = 35.0;                          % Maximal frequency in Hz
aWav.aName = 'cwt.dat';                     % File with wavelet coefficients
aWav.aRidge = 'cwtridge.dat';               % Parameter of the wavelet, time-frequency resolution (f/df)
aWav.aWavelet = 'morlet --points=1024';     % Name of the wavelet. In v.1.4 only two possibilities: 'morlet' or 'cauchy'
aWav.aPar = 2.0;                            % Parameter of the wavelet, time-frequency resolution (f/df)
% Wavelet spectrum
[aTime,aSignal,aParSig] = gwlSignalRead(1,aSource.Name,'func',['--format=ASCII --smplfreq=100 --tmin=' num2str(aSource.Tmin) ' --tmax=' num2str(aSource.Tmax) ' --to2p'],aSource.Bin);
gwlCreateAxis(aWav.aFCount,aWav.aFMin,aWav.aFMax,'lin',aWav.aFName,'Frequency');
gwlCwt(1, aSource.Bin, aWav.aFName, 2, aWav.aWavelet, aWav.aPar, aWav.aName);
[aTimeWT,aFreqWT,aCwt] = gwlConvert('3','',aWav.aName);
% Sceleton calculation
gwlExec('gwlCwtMaxLine', ['--infile=' aWav.aName ' --outfile=' aWav.aRidge ' --type=amaxt'])
fid = fopen(aWav.aRidge,'r');  
  [aTimeR,aRidge] = gwlReadSignal(fid);
fclose(fid);  
% Fourier spectrum and averaged wavelet spectrum
[aFreq1,aFour1] = gwlFourTrans('mat',aSignal,1,aParSig.aSample);
aCwtFour = sum(aCwt,2)/(2*pi);
% Plots
gwlPlotFunction(aTime,aSignal,0.07,0.75,0.9,0.20,min(aTime),max(aTime),min(aSignal),max(aSignal),'',gwlGetNotation('VSIG','T','z'),'(a) Geophone: 1, Cannel: Z');
gwlPlotImage(aTimeWT,aFreqWT, aCwt, 0.07,0.33,0.9,0.40, gwlGetNotation('TIME'), gwlGetNotation('FREQ'), ['(b) ' gwlGetNotation('VSIG','WABS','z')]);
     line(aTimeR,aRidge,'Color',gwlGetColor(0));
gwlPlotFunction(aFreq1,abs(aFour1),0.32,0.07,0.39,0.20,aWav.aFMin,aWav.aFMax,0,max(abs(aFour1)),gwlGetNotation('FREQ'),gwlGetNotation('VSIG','FM','z'),'(c)');
    hold on;   plot(aFreqWT,aCwtFour,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',2);    hold off;
delete(aSource.Bin); delete(aWav.aFName); delete(aWav.aName); delete(aWav.aRidge);
    
%---------------------------------------------------------------------------
pause(0.00001);
print -f2 -r600 -depsc VibroBuildingFig2;
print -f3 -r600 -depsc VibroBuildingFig3;
print -f4 -r600 -depsc VibroBuildingFig4;
print -f5 -r600 -depsc VibroBuildingFig5;


