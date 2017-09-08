function TwoNoised()
% TwoNoised(): In order to demonstrate the noise stability, we perform the ‘‘frequency-
% velocity’’ analysis for the example presented in TwoModes.m but add random noise to 
% the source signal. The noise level is about 5 percent related to the peak-to-peak 
% amplitude of the original signal. These noisy seismograms are plotted in Figure 1a; 
% Figures 1b,c show the correlation images.
% 
% [1] M. Kulesh, M. Holschneider, M. Ohrnberger, E. Lueck. Modeling of wave dispersion using 
%     continuous wavelet transforms II: wavelet based frequency-velocity analysis // Pure and 
%     Applied Geophysics. Vol. 165. (2008, in press). DOI 10.1007/s00024-008-0299-7
% 
% FIGURE 1. Frequency-velocity analysis of noised seismic arrivals with two interfering 
% pulses.
% 
% IMPORTANT! Before run this example, go to ./BinData and run 
% TwoNoised.bat from there. The calculation take about 3 hours.

%---------------------------------------------------------------------------
path(path, '../../mshell');

%---------------------------------------------------------------------------
figure(1);
Par.aGrad = 1;
Par.SeisName = 'BinData/TwoNoisedSig.dat';
Par.FKName1 = 'BinData/TwoNoisedArg.dat';
Par.FKName2 = 'BinData/TwoNoisedCphase.dat';
Par.Mode1 = 'BinData/TwoNoisedMod1.dat';
Par.Mode2 = 'BinData/TwoNoisedMod2.dat';
plotsData(Par);

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc TwoNoised;

%---------------------------------------------------------------------------
% Local function
%---------------------------------------------------------------------------
function plotsData(aPar)
aYmax = 50;
aVelCol = gwlGetColor(0);
fid = fopen(aPar.SeisName,'r'); [aTime,aSignal]=gwlReadSignal(fid); fclose(fid);
fid = fopen(aPar.FKName1,'r'); [aVel,aFreq,aFK1]=gwlReadSpectrum(fid); fclose(fid);
fid = fopen(aPar.FKName2,'r'); [aVel,aFreq,aFK2]=gwlReadSpectrum(fid); fclose(fid);
fid = fopen(aPar.Mode1,'r'); [aFreqM,aMode1]=gwlReadDispModel(fid); fclose(fid);
fid = fopen(aPar.Mode2,'r'); [aFreqM,aMode2]=gwlReadDispModel(fid); fclose(fid);
gwlPlotSeis(aTime,aSignal,0.07,0.78,0.9,0.2,min(aTime),max(aTime),aYmax,gwlGetNotation('TIME'),gwlGetNotation('TRN'),0,'(a)');
gwlPlotImage(aFreq,aVel,(aFK1'.^(aPar.aGrad)),0.07,0.07,0.405,0.65,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CP','F'),'');
    line(aFreqM,aMode1(:,3),'Color',aVelCol,'LineWidth',2);
    line(aFreqM,aMode2(:,3),'Color',aVelCol,'LineWidth',2);
    gwlText(230,1445,'(b)');
gwlPlotImage(aFreq,aVel,(aFK2'.^(aPar.aGrad)),0.48,0.07,0.49,0.65,gwlGetNotation('FREQ'),'','');
    line(aFreqM,aMode1(:,3),'Color',aVelCol,'LineWidth',2);
    line(aFreqM,aMode2(:,3),'Color',aVelCol,'LineWidth',2);
    gwlText(230,1445,'(c)');
    colorbar;

    
    