function RayleightZoom()
% RayleightZoom(): Two component polarization analysis usind simulated seismogram with
% Rayleight-wave arrivals.
%
% The seismograms in this example are computed for an impulsive force acting on the 
% surface of a single layer over a half-space model using the modal-summation method 
% of Herrmann (Herrmann, R. B., 1996, Computer programs in seismology: Saint Louis
% University, Version 3.0.). The large-amplitude wave package arriving is produced 
% predominantly by the fundamental-mode Ray-leigh-wave component, while the body-wave 
% contributions arrive at earlier times.
% 
% [1] M.S.Diallo, M.Holschneider, M.Kulesh, F.Scherbaum and F.Adler.
%     Characterization of seismic waves polarization attributes using continuous wavelet 
%     transforms // Preprint Series DFG SPP 1114, University of Bremen. Preprint 38 (2003).
% [2] M.S.Diallo, M.Kulesh, M.Holschneider, F.Scherbaum, F.Adler.
%     Characterization of polarization attributes of seismic waves using continuous wavelet 
%     transforms // Geophysics. V. 71. No. 3. P. V67-V77 (2006).
% [3] M.A.Kulesh, M.S.Diallo and M.Holschneider. Wavelet analysis of ellipticity, dispersion, 
%     and dissipation properties of Rayleigh waves // Acoustical Physics. V. 51. No. 4. 
%     P. 421-434 (2005).
% [4] M.S.Diallo, M.Holschneider, F.Scherbaum and M.Kulesh. Characterization of dispersive Rayleigh 
%     waves using wavelet transform // Eos. Trans. AGU, 84(46), Fall Meet. Suppl., Abstract 
%     S22B-0442 (2003).
% 
% FIGURE 1. (a) Components of the synthetic 2-C seismograms. (b) Its progressive and regressive
% wavelet transform. (c) Hodogram showing the particle motion over the entire time window.
%     
% FIGURE 2. (a) The semimajor and semiminor axes in time-frequency domain obtained from the
% present study for the synthetic seismograms in Figure 1. (b) The corresponding reciprocal 
% ellipticity. (c) The phase difference between the two components. (d) The rise angle. The 
% dashed lines correspond to the CTA method (Rene et. all, 1986, Multicomponent seismic studies 
% using complex trace analysis: Geophysics, 51, 1235-1251.), while the solid lines show the 
% values obtained from the present study.    
%
% FIGURE 3. (a) Variation of the signed ellipticity versus frequency and (b) ellipticity versus 
% time. The dashed lines show the theoretical values, while the solid lines show the values 
% obtained from the present study.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';
aRidgeName1 = 'ridge1.dat';
aRidgeName2 = 'ridge2.dat';
aElliparName = 'elli.dat';
aTmin = 0.0;
aTmax = 10;
aYmax = 500;
aImax = 2002;

%---------------------------------------------------------------------------
[aTime,aSignal,aParSig] = gwlSignalRead(2,'RayleightZoom.asc','func','--format=ASCII --chan=0,1 --istime',aSignalName,'Synthetic signal');
aSegInd1 = 129;
aSegInd2 = 1345;

gwlCreateAxis(128,0.01,4,'lin --sign=full',aFreqName,'Frequency');
gwlCwt(2, aSignalName, aFreqName, 1, 'cauchy', 5, aSpectrName);
[aTime,aFreq,aCwt] = gwlConvert('3','',aSpectrName);
AF1 = min(aFreq);
AF2 = max(aFreq);
[aTimeR,aRidgeT] = gwlSignalRead(1,'RayleightZoom.mlt','seis',' --istime --chan=0,1',aRidgeName1,'Time maximum line');
[aFreqR,aRidgeF] = gwlSignalRead(1,'RayleightZoom.mlf','seis',' --istime --chan=0,1',aRidgeName2,'Frequency maximum line');

gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=complex --filter=1']);
[aTimeE,aFreqE,aCwtE] = gwlConvert('1,2','',aElliparName);
[aTime1,aRmax1,aRmin1,aPhidiff1,aTilt1,aRatio1] = localCalcElliSignal(aSpectrName,1,aRidgeName1);
[aTime2,aRmax2,aRmin2,aPhidiff2,aTilt2,aRatio2] = localCalcElliSignal(aSignalName,2,'');
[aFreq3,aRmax3,aRmin3,aPhidiff3,aTilt3,aRatio3] = localCalcElliSignal(aSpectrName,3,aRidgeName2);
[aTime4,aRatio4] = textread('RayleightZoom.rat','%f %f');
[Vf,aFreq5]=pmtm(real(aSignal),1.3,length(aSignal),aParSig.aSample);
[Hf,aFreq5]=pmtm(imag(aSignal),1.3,length(aSignal),aParSig.aSample);
aRatio5 = Hf./Vf;

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
    line(aTimeR,aRidgeT(:,1),'Color',gwlGetColor(1));
    line(aTimeR,aRidgeT(:,2),'Color',gwlGetColor(1));
    gwlText(2.0,2,gwlGetNotation('CSIG','WABSP'));
    gwlText(2.0,-2,gwlGetNotation('CSIG','WABSM'));
gwlPlotFunction(real(aSignal),imag(aSignal),0.4,0.07,0.24,0.29,-aYmax,aYmax,-aYmax,aYmax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('MSIG','T','z'),'(c)');

%---------------------------------------------------------------------------
figure(2);
aSemiAxes = cat(1,aCwtE(:,:,1),aCwtE(:,:,2));
gwlPlotImage(aTime(1:aImax), aFreq, aSemiAxes(:,1:aImax), 0.07,0.76,0.9,0.22 ,'',gwlGetNotation('FREQ'),'(a)');
    set(gca,'YTick',AF1:AF2/4:AF2);
    set(gca,'YTickLabel',{0,AF2/4,3*AF2/4,AF2/2,0,AF2/4,3*AF2/4,AF2/2,AF2});
    gwlText(1.5,2.25,gwlGetNotation('EPAR','RMIN','WG'));
    gwlText(1.5,-2.25,gwlGetNotation('EPAR','RMAX','WG'));
    line([0,aTmax],[0,0],'Color','black');
gwlPlotFunction(aTime1(aSegInd1:aSegInd2),aRatio1(aSegInd1:aSegInd2),0.07,0.528,0.9,0.22,aTmin,aTmax,0,1,'',gwlGetNotation('EPAR','RATIO','T'),'(b)');
    hold on;   plot(aTime2,aRatio2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
gwlPlotFunction(aTime1(aSegInd1:aSegInd2),aPhidiff1(aSegInd1:aSegInd2),0.07,0.296,0.9,0.22,aTmin,aTmax,-pi,pi,'',gwlGetNotation('EPAR','PDIFF','T'),'(c)');
    hold on;   plot(aTime2,aPhidiff2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
gwlPlotFunction(aTime1(aSegInd1:aSegInd2),aTilt1(aSegInd1:aSegInd2),0.07,0.06,0.9,0.22,aTmin,aTmax,-pi,pi,gwlGetNotation('TIME'),gwlGetNotation('EPAR','TILT','T'),'(d)');
    hold on;   plot(aTime2,aTilt2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});

%---------------------------------------------------------------------------
figure(3);
aParName1 = ['1/' gwlGetNotation('EPAR','RATIO','F') '^2'];
aParName2 = ['sign(' gwlGetNotation('EPAR','PDIFF','T') ')/' gwlGetNotation('EPAR','RATIO','T') '^2'];
gwlPlotFunction(aFreq3,0.015./aRatio3.^2,0.07,0.57,0.9,0.4,0,AF2,0,25,gwlGetNotation('FREQ'),aParName1,'(a)');
    hold on;    plot(aFreq5,aRatio5,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);   hold off;
gwlPlotFunction(aTime(aSegInd1:aSegInd2),0.05.*sign(aTilt1(aSegInd1:aSegInd2))./aRatio1(aSegInd1:aSegInd2).^2,0.07,0.07,0.9,0.4,aTmin,aTmax,-50,50,gwlGetNotation('TIME'),aParName2,'(b)');
    hold on;    plot(aTime4,aRatio4,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);   hold off;

%---------------------------------------------------------------------------
figure(4);
gwlPlotFunction(aTime,real(aSignal),0.07,0.68,0.9,0.15,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('CSIG','T'),'(a)');
    hold on;    plot(aTime,imag(aSignal),'--black','LineWidth',1);    hold off;
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotImage(aTime(1:aImax), aFreq, aCwt(:,1:aImax), 0.07,0.43,0.9,0.24, gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(b)');
    line([0,aTmax],[0,0],'Color','black');
    line(aTimeR,aRidgeT(:,1),'Color',gwlGetColor(1));
    line(aTimeR,aRidgeT(:,2),'Color',gwlGetColor(1));
    gwlText(2.0,2,gwlGetNotation('CSIG','WABSP'));
    gwlText(2.0,-2,gwlGetNotation('CSIG','WABSM'));
gwlPlotFunction(real(aSignal),imag(aSignal),0.2,0.16,0.25,0.2,-aYmax,aYmax,-aYmax,aYmax,gwlGetNotation('CSIG','RET'),gwlGetNotation('CSIG','IMT'),'(c)');
gwlPlotFunction(aFreq3,0.015./aRatio3.^2,0.6,0.16,0.25,0.2,0,AF2,0,25,gwlGetNotation('FREQ'),aParName1,'(d)');
    hold on;    plot(aFreq5,aRatio5,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);   hold off;

%---------------------------------------------------------------------------
delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aRidgeName1); delete(aRidgeName2); delete(aElliparName);
clear all;
print -f1 -r600 -depsc RayleightZoomFig1;
print -f2 -r600 -depsc RayleightZoomFig2;
print -f3 -r600 -depsc RayleightZoomFig3;
print -f4 -r600 -depsc RayleightZoomFig4;


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
if(aType == 3)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=mlinef --filter=1 --mline=' aRidgeName]);
end;
gwlExec('gwlConvert',[' --infile=' aTmpName1 ' --outfile=' aTmpName1 ' --outtype=1 --comp=1,2,13,14,4 --nomess']);
[aTime,aRmax,aRmin,aPhidiff,aTilt,aRatio] = textread(aTmpName1,'%f %f %f %f %f %f');
delete(aTmpName1);

