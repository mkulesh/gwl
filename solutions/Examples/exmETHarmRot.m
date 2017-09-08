function exmETHarmRot()
% exmETHarmRot(): In this example we use synthetic data to illustrate how the polarization
% attributes can be recovered using the wavelet transform.
% 
% [1] M.S.Diallo, M.Holschneider, M.Kulesh, F.Scherbaum and F.Adler. Characterization 
%     of seismic waves polarization attributes using continuous wavelet transforms // 
%     Preprint Series DFG SPP 1114, University of Bremen. Preprint 38 (2003).
% [2] M.S.Diallo, M.Holschneider, M.Kulesh, F.Scherbaum and F.Adler
%     Characterization of the Rayleigh wave polarization attributes with continuous 
%     wavelet transform // Geophysical Research Abstracts, Vol. 5, 11237 (2003).
%
% FIGURE 1. (a) Component of the complex signal C(t), (b) its wavelet transform 
% obtained  with a regressive  Cauchy wavelet  (bottom) and  progressive  
% Cauchy wavelet (top) and (c) plot of x(t) vs y(t) components
% that displays the rotating ellipses
%
% FIGURE 2. (a) The  wavelet transform  of the  major  half-axis R(b,a) and minor half-axis r(b,a), 
% (b) The reciprocal ellipticity from the present study where averaging over the frequencies 
% was performed (solid lines) and the reciprocal ellipticity from Rene CTA (dashed line).
% (c) The phase difference and (d) the tilt angle in comparsion with Rene CTA.
%
% FIGURE 3. The  major  half-axis R(t), the minor half-axis r(t), the phase difference 
% and the tilt angle as functions of time obtained from the present 
% study (solid lines) in comparsion with Rene CTA (dashed line).

%---------------------------------------------------------------------------
path(path, '../../mshell');

aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';
aElliparName = 'elli.dat';
aYmax = 10.0;
aTmin = 0.0;
aTmax = 5.0;
aImax = 1001;

%---------------------------------------------------------------------------
gwlCreateAxis(1024,0,5.115,'lin',aTimeName,'Time');
[aTime,aSignal,aParSig] = gwlSignalGen(2,aTimeName,'harmrot','2.0,7.0,1.0,2.0,5.0,5.0,0.318',aSignalName,'Rotated harmonic function');

gwlCreateAxis(128,0.0001,15,'lin --sign=full',aFreqName,'Frequency');
gwlCwt(2, aSignalName, aFreqName, 0, 'cauchy', 4, aSpectrName);
[aTime,aFreq,aCwt] = gwlConvert('3','',aSpectrName);

gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=complex']);
[aTimeE,aFreqE,aCwtE] = gwlConvert('1,2','',aElliparName);
[aTime1,aRmax1,aRmin1,aPhidiff1,aTilt1,aRatio1,aWx1,aWy1] = calcElliSignal(aSignalName,1);
[aTime2,aRmax2,aRmin2,aPhidiff2,aTilt2,aRatio2,aWx2,aWy2] = calcElliSignal(aSignalName,2);
[aTime3,aRmax3,aRmin3,aPhidiff3,aTilt3,aRatio3,aWx3,aWy3] = calcElliSignal(aElliparName,3);

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aTime,real(aSignal),0.07,0.825,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','x'),'(a)');
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotFunction(aTime,imag(aSignal),0.07,0.68,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','z'),'');
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotImage(aTime(1:aImax), aFreq, aCwt(:,1:aImax), 0.07,0.43,0.9,0.24, gwlGetNotation('TIME'), gwlGetNotation('FREQ'), '(b)');
    line([0,aTmax],[0,0],'Color','black');
    gwlText(1.0,10,gwlGetNotation('CSIG','WABSP'));
    gwlText(1.0,-7,gwlGetNotation('CSIG','WABSM'));
gwlPlotFunction(real(aSignal),imag(aSignal),0.4,0.07,0.24,0.29,-aYmax,aYmax,-aYmax,aYmax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('MSIG','T','z'),'(c)');

%---------------------------------------------------------------------------
figure(2);
aSemiAxes = cat(1,aCwtE(:,:,1).^2,aCwtE(:,:,2).^2);
gwlPlotImage(aTime(1:aImax), aFreq, aSemiAxes(:,1:aImax), 0.07,0.76,0.9,0.22 ,'',gwlGetNotation('FREQ'),'(a)');
    AF1 = min(aFreq);
    AF2 = max(aFreq);
    set(gca,'YTick',AF1:AF2/3:AF2);
    set(gca,'YTickLabel',{0,AF2/3,2*AF2/3,0,AF2/3,2*AF2/3,AF2});
    gwlText(1.1,8,gwlGetNotation('EPAR','RMIN','WG'));
    gwlText(1.1,-4,gwlGetNotation('EPAR','RMAX','WG'));
    line([0,aTmax],[0,0],'Color','black');
gwlPlotFunction(aTime3,aRatio3,0.07,0.528,0.9,0.22,aTmin,aTmax,0,1,'',gwlGetNotation('EPAR','RATIO','T'),'(b)');
    hold on;   plot(aTime2,aRatio2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
gwlPlotFunction(aTime3,aPhidiff3,0.07,0.296,0.9,0.22,aTmin,aTmax,-pi,pi,'',gwlGetNotation('EPAR','PDIFF','T'),'(c)');
    hold on;   plot(aTime2,aPhidiff2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
gwlPlotFunction(aTime3,aTilt3,0.07,0.06,0.9,0.22,aTmin,aTmax,-pi,pi,gwlGetNotation('TIME'),gwlGetNotation('EPAR','TILT','T'),'(d)');
    hold on;   plot(aTime2,aTilt2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});

%---------------------------------------------------------------------------
figure(3);
gwlPlotFunction(aTime,real(aSignal),0.07,0.78,0.9,0.175,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('CSIG','T'),'(a)');
    hold on;   plot(aTime,imag(aSignal),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    legend(gwlGetNotation('CSIG','RET'), gwlGetNotation('CSIG','IMT'))
gwlPlotFunction(aTime1,aPhidiff1,0.07,0.60,0.9,0.165,aTmin,aTmax,-pi,pi,'',gwlGetNotation('EPAR','PDIFF','T'),'(b)');
    hold on;   plot(aTime2,aPhidiff2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
gwlPlotFunction(aTime1,aRmax1,0.07,0.42,0.9,0.165,aTmin,aTmax,0,aYmax,'',[gwlGetNotation('EPAR','RMAX','T') ', ' gwlGetNotation('EPAR','RMIN','T')],'(c)');
    hold on;   plot(aTime2,aRmax2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    hold on;   plot(aTime1,aRmin1,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime2,aRmin2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    gwlText(2.5,7.1,gwlGetNotation('EPAR','RMAX','T'));
    gwlText(2.5,3.6,gwlGetNotation('EPAR','RMIN','T'));
gwlPlotFunction(aTime1,aTilt1,0.07,0.24,0.9,0.165,aTmin,aTmax,-pi,pi,'',gwlGetNotation('EPAR','TILT','T'),'(d)');
    hold on;   plot(aTime2,aTilt2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
gwlPlotFunction(aTime1,aWx1,0.07,0.06,0.9,0.165,aTmin,aTmax,-6,6,gwlGetNotation('TIME'),gwlGetNotation('EPAR','INST','T'),'(e)');
    hold on;   plot(aTime1,aWy1,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    gwlText(2.5,4.8,gwlGetNotation('EPAR','INST','T',1));
    gwlText(2.5,-0.5,gwlGetNotation('EPAR','INST','T',2));

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aElliparName); 
clear all;

print -f1 -r600 -depsc exmETHarmRotFig1;
print -f2 -r600 -depsc exmETHarmRotFig2;
print -f3 -r600 -depsc exmETHarmRotFig3;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function [aTime,aRmax,aRmin,aPhidiff,aTilt,aRatio,aWx,aWy] = calcElliSignal(aSourceName,aType)
aTmpName1 = tempname;
if(aType == 1)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=complex']);
end;
if(aType == 2)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=rene']);
end;
if(aType == 3)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --freq=0.1,15']);
end;
gwlExec('gwlConvert',[' --infile=' aTmpName1 ' --outfile=' aTmpName1 ' --outtype=1 --comp=1,2,13,14,4,7,8 --nomess']);
[aTime,aRmax,aRmin,aPhidiff,aTilt,aRatio,aWx,aWy] = textread(aTmpName1,'%f %f %f %f %f %f %f %f');
delete(aTmpName1);

