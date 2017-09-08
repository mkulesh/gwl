function WorthBasinSynt()
% In this example we separate polarized events from a 2-C elastic synthetic example from 
% the Fort Worth Basin, U. S. A. The identification of polarized events in these synthetics 
% is complicated by the presence of considerable interfering events, such as P- and S-waves 
% with backscattered head waves. The synthetics also contain converted as well as surface
% waves. Details about the construction of the model and synthetic generation are available 
% in (Gherasim, M., C. Hoelting, and K. Marfurt, 2004, Fort Worth Basin - 2-Delastic synthetic 
% depth model: 74th Annual International Meeting, SEG, Expanded Abstracts, 1893-1896).
% 
% [1] M.S.Diallo, M.Kulesh, M.Holschneider, F.Scherbaum, F.Adler.
%     Characterization of polarization attributes of seismic waves using continuous wavelet 
%     transforms // Geophysics. V. 71. No. 3. P. V67-V77 (2006).
% [2] Diallo M.S., Kulesh M., Holschneider M. and Scherbaum F. Instantaneous polarization 
%     attributes in the time-frequency domain: application to wave field separation // Eos 
%     Trans. AGU, 85(47), Fall Meet. Suppl., Abstract S31B-1063 (2004).
% [3] M.Holschneider, M.S.Diallo and K. Kurennaya. Elliptic Properties of Surface Waves in 
%     Wavelet Domain // Book of abstracts of XXXIII Summer School - Conference "Advanced 
%     Problems in Mechanics" (June 28 - July 5, 2005, St.Petersburg, Russia). P.47.
% 
% FIGURE 1. 2-C elastic synthetic seismogram from the Fort Worth Basin. (a) The horizontal
% component. (b) The vertical component. (c) Computed reciprocal ellipticity from the CTA. 
% (d) Computed rise angle form the CTA. From last two plots we can identify the polarization 
% characteristics of the prominent wave modes present in the synthetics.
% 
% FIGURE 2 and FIGURE 3. Separation of wave modes present in the 2-C elastic synthetic seismogram 
% from the Fort Worth Basin, denoted by LH, LV, EH, and EH, with ratio = 0.15 and tilt angle = 0.7.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aSeisMax = 1000;

%---------------------------------------------------------------------------
figure(1);
fid = fopen('BinData/WorthBasinSyntsig.dat','r');    [aTime,aSignal]=gwlReadSignal(fid);  fclose(fid);
fid = fopen('BinData/WorthBasinSyntratio.dat','r');  [aTime,aRatio,aParRatio]=gwlReadSignal(fid);  fclose(fid);
fid = fopen('BinData/WorthBasinSynttilt.dat','r');   [aTime,aTilt,aParTilt]=gwlReadSignal(fid);  fclose(fid);

gwlPlotSeisV(aTime,real(aSignal),0.07,0.53,0.37,0.45,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(a)');
    localLabelHoriz;
gwlPlotSeisV(aTime,imag(aSignal),0.53,0.53,0.37,0.45,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(b)');
    localLabelVert;
gwlPlotImage(1:aParRatio.aChanCount,aTime,aRatio,0.07,0.07,0.43,0.43,gwlGetNotation('MSIG','T','x'),gwlGetNotation('TIME'),['(c) ' gwlGetNotation('EPAR','RATIO','T')]);
    colorbar;
    set(gca,'YDir','reverse');
    obj2 = get(gca,'Children');  set(obj2(1),'Position',[2.6 0.05 0]);
    set(gwlColorBar,'YTick',0:0.2:1);
    set(gwlColorBar,'YTickLabel',{0,0.2,0.4,0.6,0.8,1});
    localLabelHoriz;
    localLabelVert;
gwlPlotImage(1:aParTilt.aChanCount,aTime,aTilt,0.53,0.07,0.43,0.43,gwlGetNotation('MSIG','T','z'),'',['(d) ' gwlGetNotation('EPAR','TILT','T')]);
    colorbar;
    set(gca,'YDir','reverse');
    obj2 = get(gca,'Children');  set(obj2(1),'Position',[2.6 0.05 0]);
    set(gwlColorBar,'YTick',-1.5:0.75:1.5);
    set(gwlColorBar,'YTickLabel',{'-pi/2','-pi/4','0','pi/4','pi/2'});
    localLabelHoriz;
    localLabelVert;
colormap hsv;
 
%---------------------------------------------------------------------------
figure(2);
fid = fopen('BinData/WorthBasinSyntfilt(1).dat','r'); [aTime,aLinHor]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/WorthBasinSyntfilt(2).dat','r'); [aTime,aLinVer]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/WorthBasinSyntfilt(3).dat','r'); [aTime,aEllHor]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/WorthBasinSyntfilt(4).dat','r'); [aTime,aEllVer]=gwlReadSignal(fid); fclose(fid);

gwlPlotSeisV(aTime,real(aLinHor),0.25,0.78,0.25,0.21,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(a)');
gwlPlotSeisV(aTime,imag(aLinHor),0.53,0.78,0.25,0.21,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(b)');
gwlPlotSeisV(aTime,real(aLinVer),0.25,0.54,0.25,0.21,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(a)');
gwlPlotSeisV(aTime,imag(aLinVer),0.53,0.54,0.25,0.21,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(b)');
gwlPlotSeisV(aTime,real(aEllHor),0.25,0.30,0.25,0.21,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(a)');
gwlPlotSeisV(aTime,imag(aEllHor),0.53,0.30,0.25,0.21,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(b)');
gwlPlotSeisV(aTime,real(aEllVer),0.25,0.06,0.25,0.21,min(aTime),max(aTime),aSeisMax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('TIME'),0,'(a)');
gwlPlotSeisV(aTime,imag(aEllVer),0.53,0.06,0.25,0.21,min(aTime),max(aTime),aSeisMax,gwlGetNotation('MSIG','T','z'),gwlGetNotation('TIME'),0,'(b)');
    
%---------------------------------------------------------------------------
figure(3);
gwlPlotSeisV(aTime,real(aEllHor),0.09,0.51,0.41,0.41,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(a)');
    localLabelHoriz;
gwlPlotSeisV(aTime,imag(aEllHor),0.53,0.51,0.41,0.41,min(aTime),max(aTime),aSeisMax,' ',' ',0,'(b)');
    localLabelVert;
gwlPlotSeisV(aTime,real(aEllVer),0.09,0.07,0.41,0.41,min(aTime),max(aTime),aSeisMax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('TIME'),0,'(c)');
    localLabelHoriz;
gwlPlotSeisV(aTime,imag(aEllVer),0.53,0.07,0.41,0.41,min(aTime),max(aTime),aSeisMax,gwlGetNotation('MSIG','T','z'),' ',0,'(d)');
    localLabelVert;

%---------------------------------------------------------------------------
figure(4);
gwlPlotSeisV(aTime,real(aLinHor),0.09,0.51,0.41,0.41,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(a)');
    localLabelHoriz;
gwlPlotSeisV(aTime,imag(aLinHor),0.53,0.51,0.41,0.41,min(aTime),max(aTime),aSeisMax,' ',' ',0,'(b)');
    localLabelVert;
gwlPlotSeisV(aTime,real(aLinVer),0.09,0.07,0.41,0.41,min(aTime),max(aTime),aSeisMax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('TIME'),0,'(c)');
    localLabelHoriz;
gwlPlotSeisV(aTime,imag(aLinVer),0.53,0.07,0.41,0.41,min(aTime),max(aTime),aSeisMax,gwlGetNotation('MSIG','T','z'),' ',0,'(d)');
    localLabelVert;

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc WorthBasinSyntFig1;
print -f2 -r600 -depsc WorthBasinSyntFig2;
print -f3 -r600 -depsc WorthBasinSyntFig3;
print -f4 -r600 -depsc WorthBasinSyntFig4;


%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function localLabelHoriz()
    gwlText(143,0.7,'A');
    gwlText(74,0.25,'B');
    gwlText(92,0.53,'C');
    gwlText(146,1.8,'D');
    gwlText(60,0.73,'E');
    gwlText(65,1.08,'F');
    gwlText(90,1.37,'G');
    gwlText(116,1.75,'H');
    gwlText(55,1.85,'I');
    gwlText(113,2.08,'J');
    gwlText(151,2.68,'K');
    gwlText(34,2.80,'L');
    gwlText(43,3.46,'M');

function localLabelVert()
    gwlText(76,0.11,'N');
    gwlText(150,0.73,'O');
    gwlText(80,0.75,'P');
    gwlText(80,1.26,'Q');
    gwlText(80,1.7,'R');
    gwlText(141,1.65,'S');
    gwlText(137,2.02,'T');
