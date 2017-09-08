function RayleightFull()
% RayleightFull(): Two component polarization analysis usind simulated seismogram with
% Rayleight-wave arrivals.
%
% The seismograms in this example are computed for an impulsive force acting on the 
% surface of a single layer over a half-space model using the modal-summation method 
% of Herrmann (Herrmann, R. B., 1996, Computer programs in seismology: Saint Louis
% University, Version 3.0.). The large-amplitude wave package arriving after approximately
% 3.5 s is produced predominantly by the fundamental-mode Ray-leigh-wave component, while 
% the body-wave contributions arrive at earlier times.
% 
% [1] M.S.Diallo, M.Holschneider, M.Kulesh, F.Scherbaum and F.Adler.
%     Characterization of seismic waves polarization attributes using continuous wavelet 
%     transforms // Preprint Series DFG SPP 1114, University of Bremen. Preprint 38 (2003).
% [2] M.S.Diallo, M.Kulesh, M.Holschneider, F.Scherbaum, F.Adler.
%     Characterization of polarization attributes of seismic waves using continuous wavelet 
%     transforms // Geophysics. V. 71. No. 3. P. V67-V77 (2006).
% [3] M.Holschneider, M.S.Diallo and K. Kurennaya. Elliptic Properties of Surface Waves in 
%     Wavelet Domain // Book of abstracts of XXXIII Summer School - Conference "Advanced 
%     Problems in Mechanics" (June 28 - July 5, 2005, St.Petersburg, Russia). P.47.
% 
% FIGURE 1. (a) Components of the synthetic 2-C seismograms. (b) Its progressive and regressive
% wavelet transform. (c) Hodogram showing the particle motion over the entire time window.
% (d) Hodogram for the time window of the Rayleigh-wave arrival between 3.8 s and 8 s.
%     
% FIGURE 2. (a) The semimajor and semiminor axes in time-frequency domain obtained from the
% present study for the synthetic seismograms in Figure 1. (b) The corresponding reciprocal 
% ellipticity. (c) The phase difference between the two components. (d) The rise angle. The 
% dashed lines correspond to the CTA method (Rene et. all, 1986, Multicomponent seismic studies 
% using complex trace analysis: Geophysics, 51, 1235-1251.), while the solid lines show the 
% values obtained from the present study.    

%---------------------------------------------------------------------------
path(path, '../../mshell');
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';
aRidgeName = 'ridge.dat';
aElliparName = 'elli.dat';
aTmin = 0.0;
aTmax = 10;
aYmax = 1;
aImax = 2002;

%---------------------------------------------------------------------------
[aTime,aSignal,aParSig] = gwlSignalRead(2,'RayleightFull.sta','func','--format=STA --geo=2 --chan=0,1 --tmin=0 --tmax=10.235',aSignalName,'Synthetic signal');
aSignal = aSignal*1800;
aSegInd1 = floor(4.0*aParSig.aSample);
aSegInd2 = floor(aSegInd1+2.5*aParSig.aSample);

gwlCreateAxis(128,0.01,10,'lin --sign=full',aFreqName,'Frequency');
gwlCwt(2, aSignalName, aFreqName, 1, 'cauchy', 5, aSpectrName);
gwlExec('gwlCwtMaxLine',[' --infile=' aSpectrName ' --outfile=' aRidgeName ' --type=amaxt']);
[aTime,aFreq,aCwt] = gwlConvert('3','',aSpectrName);
fid = fopen(aRidgeName,'r');  
  [aTimeR,aRidge] = gwlReadSignal(fid);
fclose(fid);  

gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=complex --filter=1']);
[aTimeE,aFreqE,aCwtE] = gwlConvert('1,2','',aElliparName);
[aTime1,aRmax1,aRmin1,aPhidiff1,aTilt1,aRatio1] = localCalcElliSignal(aSpectrName,1,aRidgeName);
[aTime2,aRmax2,aRmin2,aPhidiff2,aTilt2,aRatio2] = localCalcElliSignal(aSignalName,2,'');

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aTime,real(aSignal),0.07,0.82,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','x'),'(a)');
      set(gca,'YTick',-aYmax:aYmax/2:aYmax);
      set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotFunction(aTime,imag(aSignal),0.07,0.68,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','z'),'');
      set(gca,'YTick',-aYmax:aYmax/2:aYmax);
      set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotImage(aTime(1:aImax), aFreq, aCwt(:,1:aImax), 0.07,0.43,0.9,0.24, gwlGetNotation('TIME'), gwlGetNotation('FREQ'), '(b)');
      line([0,aTmax],[0,0],'Color','black');
      line(aTimeR,aRidge(:,1),'Color',gwlGetColor(1));
      line(aTimeR,aRidge(:,2),'Color',gwlGetColor(1));
      gwlText(3.0,8,gwlGetNotation('CSIG','WABSP'));
      gwlText(3.0,-6,gwlGetNotation('CSIG','WABSM'));
gwlPlotFunction(real(aSignal),imag(aSignal),0.24,0.07,0.25,0.29,-aYmax,aYmax,-aYmax,aYmax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('MSIG','T','z'),'(c)');
gwlPlotFunction(real(aSignal(aSegInd1:aSegInd2)),imag(aSignal(aSegInd1:aSegInd2)),0.55,0.07,0.25,0.29,-aYmax,aYmax,-aYmax,aYmax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('MSIG','T','z'),'(d)');

%---------------------------------------------------------------------------
figure(2);
aSemiAxes = cat(1,aCwtE(:,:,1),aCwtE(:,:,2));
gwlPlotImage(aTime(1:aImax), aFreq, aSemiAxes(:,1:aImax), 0.07,0.76,0.9,0.22 ,'',gwlGetNotation('FREQ'),'(a)');
    AF1 = min(aFreq);
    AF2 = max(aFreq);
    set(gca,'YTick',AF1:AF2/4:AF2);
    set(gca,'YTickLabel',{0,AF2/4,AF2/2,3*AF2/4,0,AF2/4,AF2/2,3*AF2/4,AF2});
    gwlText(2.5,8,gwlGetNotation('EPAR','RMIN','WG'));
    gwlText(2.5,-3,gwlGetNotation('EPAR','RMAX','WG'));
    line([0,aTmax],[0,0],'Color','black');
gwlPlotFunction(aTime1,aRatio1,0.07,0.528,0.9,0.22,aTmin,aTmax,0,1,'',gwlGetNotation('EPAR','RATIO','T'),'(b)');
    hold on;   plot(aTime2,aRatio2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    legend('This study','Rene et al. (1986)');
gwlPlotFunction(aTime1,aPhidiff1,0.07,0.296,0.9,0.22,aTmin,aTmax,-pi,pi,'',gwlGetNotation('EPAR','PDIFF','T'),'(c)');
    hold on;   plot(aTime2,aPhidiff2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
    legend('This study','Rene et al. (1986)');
gwlPlotFunction(aTime1,aTilt1,0.07,0.06,0.9,0.22,aTmin,aTmax,-pi,pi,gwlGetNotation('TIME'),gwlGetNotation('EPAR','TILT','T'),'(d)');
    hold on;   plot(aTime2,aTilt2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
    legend('This study','Rene et al. (1986)');

%---------------------------------------------------------------------------
pause(0.00001);
delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aRidgeName); delete(aElliparName);
clear all;

print -f1 -r600 -depsc RayleightFullFig1;
print -f2 -r600 -depsc RayleightFullFig2;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function [aTime,aRmax,aRmin,aPhidiff,aTilt,aRatio] = localCalcElliSignal(aSourceName,aType,aRidgeName)
aTmpName1 = tempname;
if(aType == 1)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=mlinet --filter=1 --mline=' aRidgeName]);
end;
if(aType == 2)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=rene']);
end;
gwlExec('gwlConvert',[' --infile=' aTmpName1 ' --outfile=' aTmpName1 ' --outtype=1 --comp=1,2,13,14,4']);
[aTime,aRmax,aRmin,aPhidiff,aTilt,aRatio] = textread(aTmpName1,'%f %f %f %f %f %f');
delete(aTmpName1);

