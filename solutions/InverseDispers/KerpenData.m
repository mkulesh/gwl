function KerpenData()
% The experimental data (Fig. 1) consist of a 2-D shallow seismic survey (stations 
% along a line) at Kerpen, a particular site in the Lower Rhine embayment where 
% the buried scarp of a historically active fault is presumed. Several profiles of 
% 48 channels with 2-m interreceiver spacing were collected using hammer blows as 
% seismic source. From the analysis of the first arrival data, a 2-D tomographic 
% model of the investigated site was derived, however, no shear-wave velocity 
% information is available for the shallow structure. We selected a seismogram 
% profile (Fig. 1) with prominent low-frequency, high-amplitude arrivals that 
% correspond to the surface wave arrivals we intend to characterize.
% 
% We selected two subsections from the seismograms for our analysis. These 
% subsections are labelled "subsection A" and "subsection B" in Fig. 1. In 
% subsection A, the surface wave arrival consists of one coherent energy arrival 
% while in subsection B two distinct coherent energy arrivals can be observed. The 
% occurrence of these two different energy patches in the wavelet image are not 
% indicative of the typical modes of the observed Rayleigh wave as in such a case 
% a certain frequency overlap between the different modes is to be expected. In 
% the presented example the two coherent energy arrivals do not show any overlap 
% in frequency content. Therefore, we speculate that existence of these two 
% coherent arrivals in the frequency band between 5 and 40 Hz in "subsection B" is 
% most likely due to the strong attenuation characteristic of the medium under 
% investigation for a certain frequency range (around 30 Hz), which happens to lie 
% within the frequency band of the observed surface
% 
% [1] M.Holschneider, M.S.Diallo, M.Kulesh, F.Scherbaum, M.Ohrnberger and E.Lueck 
%     Characterization of dispersive surface wave using continuous wavelet transforms 
%     // Preprint Series DFG SPP 1114, University of Bremen. Preprint 56 (2004).
% 
% [2] M. Holschneider, M. S. Diallo, M. Kulesh, M. Ohrnberger, E. Lueck, F. Scherbaum 
%     Characterization of dispersive surface waves using continuous wavelet 
%     transforms // Geophysical Journal International. V. 163. No. 2. P. 463-478 (2005).
% 
% [3] M.Kulesh, M.Holschneider, M.S.Diallo, F.Scherbaum, M.Ohrnberger, E.Lueck 
%     Estimating attenuation, phase and group velocity of surface waves observed on 2D 
%     shallow seismic line using continuous wavelet transform // Proceedings of the 
%     XXXII International Summer School "Advanced Problems in Mechanics" (Russia, St. 
%     Petersburg, June 24-July 1, 2004). P. 257-262 (2004).
% 
% [4] M.Kulesh, M.Holschneider, M.S.Diallo, F.Scherbaum and M.Ohrnberger 
%     Estimating attenuation, phase and group velocity of surface waves observed on a 
%     2D shallow seismic line using continuous wavelet transform // Book of abstracts 
%     of XXXII International Summer School - Conference "Advanced Problems in 
%     Mechanics" (June 24 - July 1, 2004, St.Petersburg, Russia). P. 65-66.
% 
% [5] M.S.Diallo, M.Holschneider, M.Kulesh, F.Scherbaum, M.Ohrnberger and E.Lueck 
%     Toward improved methods of estimating attenuation, phase and group velocity of 
%     surface waves observed on shallow seismic records // Eos Trans. AGU, 85(17), Jt. 
%     Assem. Suppl., Abstract S51A-02 (2004).
% 
% [6] M.Holschneider, M.S.Diallo, M.Kulesh, F.Scherbaum and M.Ohrnberger 
%     Estimating attenuation, phase and group velocity of surface waves observed on a 
%     2D shallow seismic line using continuous wavelet transform // Geophysical 
%     Research Abstracts, Vol. 6, 03129 (2004).
% 
% FIGURE 1. Observed seismograms obtained from a shallow seismic experiment using 
% a sledgehammer as source. The distance between consecutive stations is 2 m. 
% Traces in subsections A and B are delimited by the extent of the bold lines.
% 
% FIGURE 2. Characterization of the observed surface wave as identified from the 
% traces in subsection A. (a) Comparison between the actual cross-correlations and 
% the estimated ones for selected pairs of traces. (b) The estimated phase 
% velocity (solid line) and group velocity (dashed line) and (c) the estimated 
% attenuation.
% 
% FIGURE 3. Details of the wavelet transform for signals s23(t) and s27(t) from 
% subsection B of the observed seismograms. (a) Plot of s23(t), (b) the modulus of 
% the wavelet transform of s23(t), (c) the signal s27(t) and (d) the modulus of 
% the wavelet transform of s27(t). Note that two coherent arrivals that correspond 
% to the surface wave arrivals can be observed in the frequency between 20 and 40 
% Hz.
% 
% FIGURE 4. Characterization of the observed surface wave as identified from the 
% traces in subsection B. (a) Comparison between the actual signals and the 
% estimated ones for selected traces. (b) The estimated wave number and (c) the 
% estimated attenuation.
% 
% FIGURE 5. Characterization of the low-frequency (a, b) and the high-frequency 
% (c, d) surface wave events identified in the seismograms from subsection B. (a, 
% c) The estimated phase velocity (solid line) and group velocity (dashed line) 
% and (b, d) the estimated attenuation.
% 
% FIGURE 6. (a) Comparison between the reconstructed signals (solid line) and the 
% observed seismic traces (dashed lines). (b) Comparison between the reconstructed 
% cross-correlations (solid line) and those from the observed seismograms (dashed 
% lines). In both cases a good match between the time-series can be observed. Note 
% that the cross-correlations compare much better than the signals (see text for 
% explanation).
% 
% FIGURE 7. Comparison of the slowness estimates obtained for subsection A (a) and 
% subsection B (b) from the CWT method (solid curve), Capon's highresolution f - k 
% method (contour lines) and the MUSIC algorithm (greyscaled background image).


%----------------------------------------------------------------------------
path(path, '../../mshell');

%----------------------------------------------------------------------------
figure(1);
[aTime,aSeis] = gwlSignalRead(1,'BinData/KerpenData.asc','seis',['--format=ASCII --smplfreq=4000']);
gwlPlotSeis(aTime,aSeis,0.07,0.07,0.88,0.9,min(aTime),max(aTime),2.0,gwlGetNotation('TIME'),' ',1,'');
    line([0.49,0.49],[31,41],'Color','black','LineWidth',2);
    gwlText(max(aTime)*1.01,31,'Subsec. A',90);
    line([0.49,0.49],[45,55],'Color','black','LineWidth',2);
    gwlText(max(aTime)*1.01,45,'Subsec. B',90);
    gwlText(-0.03,40,gwlGetNotation('TRN'),90);
    grid off;

%----------------------------------------------------------------------------
figure(2);
aFreqMin = 25;
aFreqMax = 45;
aAmplMax = 40;
fid = fopen('BinData/KerpenDataAoptcorr.dat','r'); [aTime,aSecAoptcorr]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/KerpenDataAcorr.dat','r'); [aTime,aSecAcorr]=gwlReadSignal(fid); fclose(fid);
gwlPlotSeis(aTime,aSecAoptcorr(:,2:4),0.07,0.5,0.88,0.4,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),'',0,'(a)');
gwlPlotSeisAdd(aTime,aSecAcorr(:,2:4),aAmplMax,0,1);
    gwlText(0.02,1.5*aAmplMax,gwlGetNotation('TSIG','T','17-18'));
    gwlText(0.02,2.5*aAmplMax,gwlGetNotation('TSIG','T','17-19'));
    gwlText(0.02,3.5*aAmplMax,gwlGetNotation('TSIG','T','17-20'));
localPlotResult('BinData/KerpenDataAoptpar.dat',2,aFreqMin,aFreqMax,0.2,'(b)','(c)');

%---------------------------------------------------------------------------
figure(3);
aChan1 = 1;
aChan2 = 5;
aAmplMax = 0.15;
fid = fopen('BinData/KerpenDataBsig.dat','r'); [aTime,aSeis]=gwlReadSignal(fid); fclose(fid);
[aTimeWT,aFreqWT,aCwtAbs] = gwlConvert('3','','BinData/KerpenDataBcwt.dat');
gwlPlotFunction(aTime,aSeis(:,aChan1), 0.07,0.87,0.88,0.12,min(aTime),max(aTime),-aAmplMax,aAmplMax,'',gwlGetNotation('SIG','T',23),'(a)');
gwlPlotImage(aTimeWT,aFreqWT,aCwtAbs(:,:,aChan1),0.07,0.615,0.88,0.25,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(b)');
    gwlText(0.3,60,gwlGetNotation('SIG','WABS',23));
gwlPlotFunction(aTime,aSeis(:,aChan2),0.07,0.40,0.88,0.12,min(aTime),max(aTime),-aAmplMax,aAmplMax,'',gwlGetNotation('SIG','T',27),'(c)');
gwlPlotImage(aTimeWT,aFreqWT,aCwtAbs(:,:,aChan2),0.07,0.145,0.88,0.25,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(d)');
    gwlText(0.4,60,gwlGetNotation('SIG','WABS',27));

%----------------------------------------------------------------------------
figure(4);
aFreqMin = 5;
aFreqMax = 100;
aAmplMax = 0.15;
fid = fopen('BinData/KerpenDataBsig.dat','r'); [aTime,aSecBsig]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/KerpenDataB1optsig.dat','r'); [aTime,aSecBopt]=gwlReadSignal(fid); fclose(fid);
gwlPlotSeis(aTime,aSecBopt,0.07,0.5,0.88,0.4,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),'',0,'(a)');
gwlPlotSeisAdd(aTime,aSecBsig,aAmplMax,0,1);
localPlotResult('BinData/KerpenDataB1optpar.dat',1,aFreqMin,aFreqMax,0.2,'(b)','(c)');

%----------------------------------------------------------------------------
figure(5);
localPlotResult('BinData/KerpenDataB2optpar.dat',2,20,30,0.57,'(a)','(b)');
localPlotResult('BinData/KerpenDataB3optpar.dat',2,30,45,0.27,'(c)','(d)');

%----------------------------------------------------------------------------
figure(6);
aAmplMax = 0.15;
gwlSignalRead(1,'BinData/KerpenData.asc','seis','--chan=23 --smplfreq=4000 --to2p','BinData/KerpenDataBsour.dat','Channel 23');
gwlExec('gwlDiffeoDisp',[' --infile=BinData/KerpenDataBsour.dat --outfile=BinData/KerpenDataB2optsig.dat --model=BinData/KerpenDataB2optpar.dat --step=4 --dist=2']);
gwlExec('gwlDiffeoDisp',[' --infile=BinData/KerpenDataBsour.dat --outfile=BinData/KerpenDataB3optsig.dat --model=BinData/KerpenDataB3optpar.dat --step=4 --dist=2']);
fid = fopen('BinData/KerpenDataBsig.dat','r'); [aTime,aSecBsig]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/KerpenDataB2optsig.dat','r'); [aTime,aSecB2optsig]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/KerpenDataB3optsig.dat','r'); [aTime,aSecB3optsig]=gwlReadSignal(fid); fclose(fid);
gwlPlotSeis(aTime,aSecB2optsig(:,2:5)+aSecB3optsig(:,2:5),0.07,0.55,0.88,0.4,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),'',0,'(a)');
gwlPlotSeisAdd(aTime,aSecBsig(:,2:5),aAmplMax,0,1);
    gwlText(0.01,1.5*aAmplMax,gwlGetNotation('SIG','T',24));
    gwlText(0.01,2.5*aAmplMax,gwlGetNotation('SIG','T',25));
    gwlText(0.01,3.5*aAmplMax,gwlGetNotation('SIG','T',26));
    gwlText(0.01,4.5*aAmplMax,gwlGetNotation('SIG','T',27));
aAmplMax = 2;
clear aTime;
fid = fopen('BinData/KerpenDataBcorr.dat','r'); [aTime,aSecBcorr]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/KerpenDataB2optcorr.dat','r'); [aTime,aSecB2optcorr]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/KerpenDataB3optcorr.dat','r'); [aTime,aSecB3optcorr]=gwlReadSignal(fid); fclose(fid);
gwlPlotSeis(aTime,aSecB2optcorr(:,2:5)+aSecB3optcorr(:,2:5),0.07,0.07,0.88,0.4,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),'',0,'(b)');
gwlPlotSeisAdd(aTime,aSecBcorr(:,2:5),aAmplMax,0,1);
    gwlText(0.02,1.5*aAmplMax,gwlGetNotation('TSIG','T','23-24'));
    gwlText(0.02,2.5*aAmplMax,gwlGetNotation('TSIG','T','23-25'));
    gwlText(0.02,3.5*aAmplMax,gwlGetNotation('TSIG','T','23-26'));
    gwlText(0.02,4.5*aAmplMax,gwlGetNotation('TSIG','T','23-27'));

%----------------------------------------------------------------------------
figure(7);
fid = fopen('BinData/KerpenDataAoptpar.dat','r');  [aFreq1,aOptPar1]=gwlReadDispModel(fid); fclose(fid);
fid = fopen('BinData/KerpenDataB2optpar.dat','r'); [aFreq2,aOptPar2]=gwlReadDispModel(fid); fclose(fid);
fid = fopen('BinData/KerpenDataB3optpar.dat','r'); [aFreq3,aOptPar3]=gwlReadDispModel(fid); fclose(fid);
localPlotMUSIC('MUSIC/music.16-19.xyz','MUSIC/capon.16-19.xyz',0.55,'(a)');
    line(aFreq1,fliplr(1./aOptPar1(:,3)),'Color','r','LineWidth',1.5);
    legend('Subsec. A');
localPlotMUSIC('MUSIC/music.23-27.xyz','MUSIC/capon.23-27.xyz',0.07,'(b)');
    line(aFreq2,fliplr(1./aOptPar2(:,3)),'Color','r','LineWidth',1.5);
    line(aFreq3,fliplr(1./aOptPar3(:,3)),'Color','r','LineWidth',1.5);
    legend('Subsec. B');

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc KerpenDataFig1;
print -f2 -r600 -depsc KerpenDataFig2;
print -f3 -r600 -depsc KerpenDataFig3;
print -f4 -r600 -depsc KerpenDataFig4;
print -f5 -r600 -depsc KerpenDataFig5;
print -f6 -r600 -depsc KerpenDataFig6;
print -f7 -r600 -depsc KerpenDataFig7;


%----------------------------------------------------------------------------
% Local function
%----------------------------------------------------------------------------
function localPlotResult(aResFile,aType,aFreqMin,aFreqMax,aButtom,aT1,aT2)
fid = fopen(aResFile,'r'); [aFreqDef,aModelDef]=gwlReadDispModel(fid); fclose(fid);
if(aType == 1)
    gwlPlotFunction(aFreqDef,aModelDef(:,1),0.07,aButtom,0.4,0.22,aFreqMin,aFreqMax,0,1.0,gwlGetNotation('FREQ'),gwlGetNotation('DISP','WN','F'),aT1);
else
    gwlPlotFunction(aFreqDef,aModelDef(:,3),0.07,aButtom,0.4,0.22,aFreqMin,aFreqMax,0,400,gwlGetNotation('FREQ'),gwlGetNotation('DISP','VEL'),aT1);
    hold on;    plot(aFreqDef,aModelDef(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend(gwlGetNotation('DISP','CP','F'),gwlGetNotation('DISP','CG','F'));
end;
gwlPlotFunction(aFreqDef,aModelDef(:,5),0.55,aButtom,0.4,0.22,aFreqMin,aFreqMax,0,0.4,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),aT2);


%----------------------------------------------------------------------------
function localPlotMUSIC(aMusic,aCapon,aVertCoor,aTitle)
indi1 = 100;
indi2 = 350;
indj1 = 85;
indj2 = 187;
CScale = 10;
YSize = 420;
[AMusic_x,AMusic_y,AMusic_z]=textread(aMusic,'%f %f %f');
[ACapon_x,ACapon_y,ACapon_z]=textread(aCapon,'%f %f %f');
FullSize = length(AMusic_x);
XSize = FullSize/YSize;
for j=1:XSize
    for i=1:YSize 
        k = (j-1)*YSize + i;
        AMXf(j) = AMusic_x(k);
        AMYf(YSize-i+1) = AMusic_y(k);
        AMZf(YSize-i+1,j) = AMusic_z(k);
        ACZf(YSize-i+1,j) = ACapon_z(k);
    end;
end;
for j=indj1:indj2
    for i=indi1:indi2
        AMX(j-indj1+1) = AMXf(j);
        AMY(i-indi1+1) = AMYf(i);
        AMZ(i-indi1+1,j-indj1+1) = AMZf(i,j);
        ACZ(i-indi1+1,j-indj1+1) = ACZf(i,j);
    end;
end;
gwlPlotSurface(AMX,fliplr(1./AMY),flipud(AMZ),2,0,50,0.07,aVertCoor,0.9,0.42,gwlGetNotation('FREQ'),['1/' gwlGetNotation('DISP','CP','F')],aTitle);
    hold on;    [C,h] = contour(AMX,fliplr(1./AMY),flipud(ACZ/CScale),5);  hold off;
