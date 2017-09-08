function exmETHarmPhase()
% exmETHarmPhase(): Two component polarization analysis usind synthetic signal
% 
% We consideras an example a sinusoidal waveform with a single frequency 
% component, constant amplitudes, and smooth phase changes. This signal is not 
% based on a physical model; our aim is to show the time-frequency resolution of 
% our method with reference to the analysis of the phase difference. We define the 
% complex synthetic signal having a small phase difference between its two 
% components as follows:
%   Z(t) = S_north(t) + i S_east(t),
%   S_north(t) = R sin(2 pi f_0 t + sin(f_1 t)),
%   S_east(t) = r sin (2 po f_0 t),
% where f_0 = 10 Hz, R = 1.5, r = 1, and f_1 = 2.
% 
% [1] M. Kulesh and M. Nose and M. Holschneider and K. Yumoto.
%     Polarization analysis of a Pi2 pulsation using continuous wavelet 
%     transform // Earth Planets Space, V. 59 (2007)
% 
% FIGURE 1. Time-frequency spectrum of the complex synthetic signal: (a) north 
% (solid line) and east (dotted line) components of the signal, (b) the modulus
% and (b) the phase of the wavelet spectrum.
% 
% FIGURE 2. Time-frequency polarization spectrum of the complex synthetic signal: 
% (a) the ellipticity ratio, (b) the tilt angle and (c) the phase difference
% between S_north(t) and S_east(t) components. 


%---------------------------------------------------------------------------
path(path, '../../mshell');

aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';
aElliparName = 'elli.dat';
aYmax = 2.0;

%---------------------------------------------------------------------------
gwlCreateAxis(1024,0,1.7,'lin',aTimeName,'Time');
[aTime,aSignal,aParSig] = gwlSignalGen(2,aTimeName,'harmphase','1.5,1,10,2',aSignalName,'Harmonic function with harmonic phase');

gwlCreateAxis(265,0.001,30,'lin --sign=full',aFreqName,'Freqency');
gwlCwt(2, aSignalName, aFreqName, 1, 'cauchy', 21, aSpectrName);
[aTime,aFreq,aCwt] = gwlConvert('3,4','--degree',aSpectrName);

gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=complex --filter=60']);
[aTimeE,aFreqE,aCwtE] = gwlConvert('4,14,13','--degree',aElliparName);

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aTime,real(aSignal),0.07,0.70,0.684,0.2,min(aTime),max(aTime),-aYmax,aYmax,'',gwlGetNotation('CSIG','T'),'(a)');
    hold on;    plot(aTime,imag(aSignal),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;

gwlPlotImage(aTime,aFreq,aCwt(:,:,1),0.07,0.49,0.80,0.2,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('CSIG','WABS')]);
   colorbar;
   grid on;

aCwt(1,1,2) = 180;  aCwt(1,2,2) = -180;
gwlPlotImage(aTime,aFreq,aCwt(:,:,2),0.07,0.28,0.80,0.2,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('CSIG','WARG')]);
   colorbar;
   grid on;

%---------------------------------------------------------------------------
figure(2);

aCwtE(1,1,1) = 0;  aCwtE(1,2,1) = 1;
gwlPlotImage(aTimeE,aFreqE,aCwtE(:,:,1),0.07,0.59,0.80,0.2,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('EPAR','RATIO','WG')]);
   colorbar;
   grid on;

aCwtE(1,1,2) = 90;  aCwtE(1,2,2) = -90;
gwlPlotImage(aTimeE,aFreqE,aCwtE(:,:,2),0.07,0.38,0.80,0.2,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','TILT','WG')]);
   colorbar;
   grid on;

aCwtE(1,1,3) = 180;  aCwtE(1,2,3) = -180;
gwlPlotImage(aTimeE,aFreqE,aCwtE(:,:,3),0.07,0.17,0.80,0.2,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('EPAR','PDIFF','WG')]);
   colorbar;
   grid on;

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aElliparName); 
clear all;
 
print -f1 -r600 -depsc exmETHarmPhaseFig1;
print -f2 -r600 -depsc exmETHarmPhaseFig2;
