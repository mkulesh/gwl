function TwoModes()
% TwoModes(): In this example, a synthetic seismogram with two modes with different
% wave numbers but with the same frequency content is analyzed. We consider the situation
% without attenuation.
% Given two different phase velocities, from which the wavenumbers are computed, a
% propagation modeling was performed; the synthetic traces are formed then by adding 
% the pulses obtained from the inverse Fourier transform. In this manner, seven 
% seismic traces were generated to simulate the observation at seven successive stations. 
% These traces are presented in Figure 1a. The phase velocities used for this seismogram 
% generation are plotted in Figures 1b,c as solid curves; fundamentally this condition 
% describes first symmetric and asymmetric modes of Lamb wave.
% 
% [1] M. Kulesh, M. Holschneider, M. Ohrnberger, E. Lueck. Modeling of wave dispersion using 
%     continuous wavelet transforms II: wavelet based frequency-velocity analysis // Pure and 
%     Applied Geophysics. Vol. 165. (2008, in press). DOI 10.1007/s00024-008-0299-7
% 
% FIGURE 1. Frequency-velocity analysis of seismic arrivals that consist of two interfering 
% pulses of different dispersion characteristics: (a) The synthetic seismogram, correlation 
% spectrum using (b) real-valued wavelet phases and (c) complex-valued phases.
% 
% FIGURE 2. Time-frequency representation of two signals used in the frequency-velocity 
% analysis in Figure 1: (a) Two analyzed signals, (b),(c) real-valued phases Wk(t,f) and 
% (d),(e) Fourier spectra of these signals.
% 
% IMPORTANT! Before run this example, go to ./BinData and run 
% TwoModes.bat from there. The calculation take about 3 hours.


%---------------------------------------------------------------------------
path(path, '../../mshell');

%---------------------------------------------------------------------------
figure(1);
Par.aGrad = 1;
Par.SeisName = 'BinData/TwoModesSig.dat';
Par.FKName1 = 'BinData/TwoModesArg.dat';
Par.FKName2 = 'BinData/TwoModesCphase.dat';
Par.Mode1 = 'BinData/TwoModesMod1.dat';
Par.Mode2 = 'BinData/TwoModesMod2.dat';
Par.Spectr1 = 'BinData/TwoModesCwtC01.dat';
Par.Spectr2 = 'BinData/TwoModesCwtC06.dat';
evalPercent(Par.SeisName,7)
plotsData(Par);

%---------------------------------------------------------------------------
figure(2);
aYmax = 50;
fid = fopen(Par.SeisName,'r'); [aTime,aSignal]=gwlReadSignal(fid); fclose(fid);
fid = fopen(Par.Spectr1,'r'); [aTime,aFreq,aCwt1]=gwlReadSpectrum(fid); fclose(fid);
fid = fopen(Par.Spectr2,'r'); [aTime,aFreq,aCwt2]=gwlReadSpectrum(fid); fclose(fid);
aSignal1(:,1) = aSignal(:,1);
aSignal1(:,2) = aSignal(:,7);
[aFreqFT,aFour] = gwlFourTrans('mat',aSignal1,2,2047);
gwlPlotSeis(aTime,aSignal1*0.7,0.07,0.74,0.9,0.25,min(aTime),max(aTime),aYmax,'',gwlGetNotation('TRN'),0,'(a)');
gwlPlotImage(aTime,aFreq,aCwt1(:,:,2),0.07,0.53,0.9,0.20,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('SIG','WARG',1)]);
gwlPlotImage(aTime,aFreq,aCwt2(:,:,2),0.07,0.32,0.9,0.20,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('SIG','WARG',2)]);
gwlPlotFunction(aFreqFT,abs(aFour(:,2)),0.07,0.07,0.4,0.2,0,250,0,max(abs(aFour(:,2))),gwlGetNotation('FREQ'),gwlGetNotation('SIG','FM',1),'(d)');
gwlPlotFunction(aFreqFT,abs(aFour(:,1)),0.57,0.07,0.4,0.2,0,250,0,max(abs(aFour(:,1))),gwlGetNotation('FREQ'),gwlGetNotation('SIG','FM',2),'(e)');

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc TwoModesFig1;
print -f2 -r600 -depsc TwoModesFig2;

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

function aVal = evalPercent(aName, aPer)
fid = fopen(aName,'r'); [aTime,aSignal]=gwlReadSignal(fid); fclose(fid);
aVal = 100*aPer/abs(max(max(aSignal)) - min(min(aSignal)));
     

    
    