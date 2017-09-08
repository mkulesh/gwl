function SynthSigB()
% In order to assess the performance of our method [2] in characterizing waves with 
% mixed coherent arrivals, we simulate seismic traces from two interfering Ricker 
% pulses with different dispersion and attenuation characteristics as shown in 
% Fig. 1. Given the phase velocity, from which the wavenumber is computed and the 
% frequency-dependent attenuation of each pulse, propagation modelling was 
% performed; the synthetic traces are then formed by adding the pulses obtained 
% from the inverse Fourier transform. In this manner, five seismic traces were 
% generated to simulate the observation at five successive stations.
% 
% From the synthesized traces at the different stations, we seek to recover the 
% phase velocity, the group velocity and the attenuation function of each pulse by 
% successive minimization of the different cost functions. For the parameter 
% vectorswe use polynomials of order 5 at most. The reason is that the parameters 
% associated with higher-order terms are negligible in comparison to those from 
% the lower terms.
% 
% [1] M.Holschneider, M.S.Diallo, M.Kulesh, F.Scherbaum, M.Ohrnberger and E.Lueck 
%     Characterization of dispersive surface wave using continuous wavelet transforms // 
%     Preprint Series DFG SPP 1114, University of Bremen. Preprint 56 (2004).
% 
% [2] M. Holschneider, M. S. Diallo, M. Kulesh, M. Ohrnberger, E. Lueck, F. 
%     Scherbaum Characterization of dispersive surface waves using continuous wavelet 
%     transforms // Geophysical Journal International. V. 163. No. 2. P. 463-478 
%     (2005).
% 
% [3] M.Kulesh, M.Holschneider, M.S.Diallo, F.Scherbaum and M.Ohrnberger 
%     Estimating attenuation, phase and group velocity of surface waves observed on a 
%     2D shallow seismic line using continuous wavelet transform // Book of abstracts 
%     of XXXII International Summer School - Conference "Advanced Problems in 
%     Mechanics" (June 24 - July 1, 2004, St.Petersburg, Russia). P. 65-66.
% 
% FIGURE 1. Modelling of a seismic arrivals that consist of two interfering Ricker 
% wavelets of different attenuation and dispersion characteristics. (a) The phase 
% velocities, (b) the group velocities, (c) the corresponding frequency-dependent 
% attenuation for each pulse, (d) overlapping wave trains propagated at different 
% stations and (e) the Fourier spectra of the different traces in (d).
% 
% FIGURE 2. Application of the wavelet optimization based on the modulus on s 1(t) 
% and s 2(t) to obtain the initial values. (a) Plot of s_1(t), (b) the modulus of 
% the wavelet transform of s_1(t). (c) the signal s_2(t) (dashed line) is compared 
% to the one obtained from the inverse wavelet transform of the estimated wavelet 
% coefficient (solid line) using derived parameters. Because the optimization is 
% based only on the modulus, the phase of the signal cannot be predicted 
% correctly. (d) The modulus of the wavelet transform of s_2(t) and (e) the 
% modulus of the wavelet transform of the estimated signal in (c).
% 
% FIGURE 3. Comparison between the original signals (dashed lines) with those 
% produced at the end of the signal optimization (solid lines). (b) Comparison 
% between the derived effective wavenumber (solid line) with those of the 
% individual pulses (dashed line). (c) Comparison between the estimated frequency-
% dependent attenuation (solid line) with those of the individual pulses (dashed 
% line). Even though a good fit of the waveforms has been achieved with this 
% optimization, the inverted parameters do not reveal explicitly the 
% characteristics of the individual pulses.
% 
% FIGURE 4. Characterization of the low-frequency pulse. (a) The modulus of the 
% wavelet transform for the cross-correlation between s_0(t) and s_3(t), (b) the 
% corresponding phase image. (c) The estimated modulus of the wavelet transform 
% for the cross-correlation between s_0(t) and s_3(t) and (d) the corresponding 
% estimated phase image. (e) Original phase velocity (dashed line) and the 
% estimated phase velocity (solid line), (f) original group velocity (dashed line) 
% and the estimated group velocity (solid line) and (g) original attenuation 
% (dashed line) and the estimated attenuation (solid line). Indeed the estimated 
% phase velocity, group velocity and the attenuation fit the initial model curves 
% very well.
% 
% FIGURE 5. Characterization of the high-frequency pulse. (a) The modulus of the 
% wavelet transform for the cross-correlation between s_0(t) and s_3(t) and (b) 
% the corresponding phase image. (c) The estimated modulus of the wavelet 
% transform for the cross-correlation between s_0(t) and s_3(t) and (d) the 
% corresponding estimated phase image. (e) Original phase velocity (dashed line) 
% and the estimated phase velocity (solid line), (f) original group velocity 
% (solid line) and the estimated group velocity (dashed line) and (g) original 
% attenuation (dashed line) and the estimated attenuation (solid line). As for the 
% low-frequency pulse, here again we observe that the estimated phase velocity, 
% group velocity and attenuation fit the initial model curves very well.
% 
% FIGURE 6. Comparison between the modelled cross-correlations (dashed lines) with 
% those obtained adding the cross-correlations derived form the inverse wavelet 
% transform (solid lines). The agreement of the time-series is very good both in 
% terms of phase and amplitude.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aLimits.aVelMin = 900;
aLimits.aVelMax = 2000;
aLimits.aAtnMax = 0.002;
fid = fopen('BinData/SynthSigBMod1.dat','r'); [aFreq,aModel1]=gwlReadDispModel(fid); fclose(fid);
fid = fopen('BinData/SynthSigBMod2.dat','r'); [aFreq,aModel2]=gwlReadDispModel(fid); fclose(fid);
fid = fopen('BinData/SynthSigBSig.dat','r'); [aTime,aSeis,aParSeis]=gwlReadSignal(fid); fclose(fid);

%---------------------------------------------------------------------------
figure(1);
aLimits.aFreqMin = 0;
aLimits.aFreqMax = 100;
aAmplMax = 1.4;
[aFreqFT,aFour] = gwlFourTrans('mat',aSeis(:,:),2,aParSeis.aSample);
gwlPlotFunction(aFreq,aModel1(:,3),0.07,0.71,0.25,0.25,aLimits.aFreqMin,aLimits.aFreqMax,aLimits.aVelMin,aLimits.aVelMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CP','F'),'(a)');
    hold on;    plot(aFreq,aModel2(:,3),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend('Mode 1','Mode 2');
gwlPlotFunction(aFreq,aModel1(:,4),0.39,0.71,0.25,0.25,aLimits.aFreqMin,aLimits.aFreqMax,aLimits.aVelMin,aLimits.aVelMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CG','F'),'(b)');
    hold on;    plot(aFreq,aModel2(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend('Mode 1','Mode 2');
gwlPlotFunction(aFreq,aModel1(:,5),0.70,0.71,0.25,0.25,aLimits.aFreqMin,aLimits.aFreqMax,0,aLimits.aAtnMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),'(c)');
    hold on;    plot(aFreq,aModel2(:,5),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend('Mode 1','Mode 2');
gwlPlotSeis(aTime,aSeis,0.07,0.34,0.88,0.31,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),'',0,'(d)');
    gwlText(0.01,1*aAmplMax,gwlGetNotation('SIG','T',0));
    gwlText(0.01,2*aAmplMax,gwlGetNotation('SIG','T',1));
    gwlText(0.01,3*aAmplMax,gwlGetNotation('SIG','T',2));
    gwlText(0.01,4*aAmplMax,gwlGetNotation('SIG','T',3));
    gwlText(0.01,5*aAmplMax,gwlGetNotation('SIG','T',4));
    grid off;
gwlPlotFunction(aFreqFT,abs(aFour),0.32,0.07,0.4,0.2,aLimits.aFreqMin,aLimits.aFreqMax,0,max(max(abs(aFour))),gwlGetNotation('FREQ'),'','(e)');
    
%---------------------------------------------------------------------------
figure(2);
aChan1 = 2;
aChan2 = 3;
[aTimeWT,aFreqWT,aCwtAbs] = gwlConvert('3','','BinData/SynthSigBCwt.dat');
[aTimeWT1,aFreqWT1,aCwtAbs1] = gwlConvert('3','','BinData/SynthSigBOpt1sp.dat');
[aTime1,aOptSig1]=gwlIwt(1, 'BinData/SynthSigBOpt1sp.dat');
gwlPlotFunction(aTime,aSeis(:,aChan1), 0.07,0.87,0.88,0.12,min(aTime),max(aTime),-aAmplMax,aAmplMax,'',gwlGetNotation('SIG','T',1),'(a)');
gwlPlotImage(aTimeWT,aFreqWT,aCwtAbs(:,:,aChan1),0.07,0.675,0.88,0.18,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('SIG','WABS',1)]);
gwlPlotFunction(aTime1,aOptSig1,0.07,0.46,0.88,0.12,min(aTime1),max(aTime1),-aAmplMax,aAmplMax,'',gwlGetNotation('SIG','T',2),'(c)');
    hold on;    plot(aTime,aSeis(:,aChan2),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
%    legend('Estimated signal','Theor. signal');
gwlPlotImage(aTimeWT,aFreqWT,aCwtAbs(:,:,aChan2),0.07,0.265,0.88,0.18,'',gwlGetNotation('FREQ'),['(d) ' gwlGetNotation('SIG','WABS',2)]);
gwlPlotImage(aTimeWT1,aFreqWT1,aCwtAbs1,0.07,0.07,0.88,0.18,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(e) ' gwlGetNotation('SIG','WABS','n')]);

%---------------------------------------------------------------------------
figure(3);
aLimits.aFreqMin = 5;
aLimits.aFreqMax = 80;
fid = fopen('BinData/SynthSigBOpt2sig.dat','r'); [aTime2,aOptSig2]=gwlReadSignal(fid); fclose(fid);
gwlPlotSeis(aTime2,aOptSig2,0.07,0.4,0.88,0.55,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),'',0,'(a)');
gwlPlotSeisAdd(aTime,aSeis,aAmplMax,0,1);
    gwlText(0.01,1*aAmplMax,gwlGetNotation('SIG','T',0));
    gwlText(0.01,2*aAmplMax,gwlGetNotation('SIG','T',1));
    gwlText(0.01,3*aAmplMax,gwlGetNotation('SIG','T',2));
    gwlText(0.01,4*aAmplMax,gwlGetNotation('SIG','T',3));
    gwlText(0.01,5*aAmplMax,gwlGetNotation('SIG','T',4));
    grid off;
localPlotResult('BinData/SynthSigBOpt2.dat',aLimits,3,1,'(b)','','(c)');

%---------------------------------------------------------------------------
figure(4);
aLimits.aFreqMin = 5;
aLimits.aFreqMax = 40;
aInd1=1;
aInd2=1023;
[aTimeWT,aFreqWT,aCwt] = gwlConvert('3,5', '--filter=5 --chan=3','BinData/SynthSigBOpt3corr.dat');
[aTimeWT1,aFreqWT1,aCwt1] = gwlConvert('3,5', '--filter=5 --chan=3','BinData/SynthSigBOpt3sp.dat');
gwlPlotImage(aTimeWT(aInd1:aInd2),aFreqWT,aCwt(:,aInd1:aInd2,1),0.07,0.64,0.44,0.24,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('TSIG','WABS','0-3')]);
gwlPlotImage(aTimeWT(aInd1:aInd2),aFreqWT,aCwt(:,aInd1:aInd2,2),0.515,0.64,0.44,0.24,'','',['(b) ' gwlGetNotation('TSIG','WARG','0-3')]);
gwlPlotImage(aTimeWT1(aInd1:aInd2),aFreqWT1,aCwt1(:,aInd1:aInd2,1),0.07,0.38,0.44,0.24,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('TSIG','WABS','n')]);
gwlPlotImage(aTimeWT1(aInd1:aInd2),aFreqWT1,aCwt1(:,aInd1:aInd2,2),0.515,0.38,0.44,0.24,gwlGetNotation('TIME'),'',['(d) ' gwlGetNotation('TSIG','WARG','n')]);
localPlotResult('BinData/SynthSigBOpt3.dat',aLimits,1,2,'(e)','(f)','(g)');

%---------------------------------------------------------------------------
figure(5);
aLimits.aFreqMin = 45;
aLimits.aFreqMax = 80;
aInd1=1;
aInd2=1023;
[aTimeWT,aFreqWT,aCwt] = gwlConvert('3,5', '--filter=5 --chan=3','BinData/SynthSigBOpt4corr.dat');
[aTimeWT1,aFreqWT1,aCwt1] = gwlConvert('3,5', '--filter=5 --chan=3','BinData/SynthSigBOpt4sp.dat');
gwlPlotImage(aTimeWT(aInd1:aInd2),aFreqWT,aCwt(:,aInd1:aInd2,1),0.07,0.64,0.44,0.24,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('TSIG','WABS','0-3')]);
gwlPlotImage(aTimeWT(aInd1:aInd2),aFreqWT,aCwt(:,aInd1:aInd2,2),0.515,0.64,0.44,0.24,'','',['(b) ' gwlGetNotation('TSIG','WARG','0-3')]);
gwlPlotImage(aTimeWT1(aInd1:aInd2),aFreqWT1,aCwt1(:,aInd1:aInd2,1),0.07,0.38,0.44,0.24,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('TSIG','WABS','n')]);
gwlPlotImage(aTimeWT1(aInd1:aInd2),aFreqWT1,aCwt1(:,aInd1:aInd2,2),0.515,0.38,0.44,0.24,gwlGetNotation('TIME'),'',['(d) ' gwlGetNotation('TSIG','WARG','n')]);
localPlotResult('BinData/SynthSigBOpt4.dat',aLimits,2,2,'(e)','(f)','(g)');

%---------------------------------------------------------------------------
figure(6);
aAmplMax = 44;
fid = fopen('BinData/SynthSigBSigCorr.dat','r'); [aTimeCorr,aCorrSig]=gwlReadSignal(fid); fclose(fid);
[aTime3,aOptSig3]=gwlIwt(1, 'BinData/SynthSigBOpt3sp.dat');
[aTime4,aOptSig4]=gwlIwt(1, 'BinData/SynthSigBOpt4sp.dat');
gwlPlotSeis(aTimeCorr,aOptSig3(:,2:5)+aOptSig4(:,2:5),0.07,0.25,0.88,0.5,min(aTimeCorr),max(aTimeCorr),aAmplMax,gwlGetNotation('TIME'),'',0,'');
gwlPlotSeisAdd(aTimeCorr,aCorrSig(:,2:5),aAmplMax,0,1);
    gwlText(0.02,1*aAmplMax,gwlGetNotation('TSIG','T','0-1'));
    gwlText(0.02,2*aAmplMax,gwlGetNotation('TSIG','T','0-2'));
    gwlText(0.02,3*aAmplMax,gwlGetNotation('TSIG','T','0-3'));
    gwlText(0.02,4*aAmplMax,gwlGetNotation('TSIG','T','0-4'));
    grid off;

%---------------------------------------------------------------------------
pause(0.00001);
clear all;

print -f1 -r600 -depsc SynthSigBFig1;
print -f2 -r600 -depsc SynthSigBFig2;
print -f3 -r600 -depsc SynthSigBFig3;
print -f4 -r600 -depsc SynthSigBFig4;
print -f5 -r600 -depsc SynthSigBFig5;
print -f6 -r600 -depsc SynthSigBFig6;

%----------------------------------------------------------------------------
% Local function
%----------------------------------------------------------------------------
function localPlotResult(aResFile,aLimits,aMode,aType,aT1,aT2,aT3)
% localPlotResult: Reading and plotting estimated velocities and attenuation
fid = fopen('BinData/SynthSigBMod1.dat','r'); [aFreq1,aModel1]=gwlReadDispModel(fid); fclose(fid);
fid = fopen('BinData/SynthSigBMod2.dat','r'); [aFreq2,aModel2]=gwlReadDispModel(fid); fclose(fid);
fid = fopen(aResFile,'r'); [aFreqDef,aModelDef]=gwlReadDispModel(fid); fclose(fid);
ASy = 0.07;
if(aType == 1)
    gwlPlotFunction(aFreqDef,aModelDef(:,1),0.23,ASy,0.25,0.22,aLimits.aFreqMin,aLimits.aFreqMax,0,0.1,gwlGetNotation('FREQ'),gwlGetNotation('DISP','WN','F'),aT1);
    if(aMode == 1 | aMode == 3)
        hold on;    plot(aFreq1,aModel1(:,1),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    end;
    if(aMode == 2 | aMode == 3)
        hold on;    plot(aFreq2,aModel2(:,1),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    end;
    AtnX1 = 0.53;
    AtnX2 = 0.25;
else
    gwlPlotFunction(aFreqDef,aModelDef(:,3),0.07,ASy,0.25,0.22,aLimits.aFreqMin,aLimits.aFreqMax,aLimits.aVelMin,aLimits.aVelMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CP','F'),aT1);
    if(aMode == 1 | aMode == 3)
        hold on;    plot(aFreq1,aModel1(:,3),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
%        legend('Estimated velocity','Theor. velocity');
    end;
    if(aMode == 2 | aMode == 3)
        hold on;    plot(aFreq2,aModel2(:,3),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
%        legend('Estimated velocity','Theor. velocity');
    end;
    gwlPlotFunction(aFreqDef,aModelDef(:,4),0.39,ASy,0.25,0.22,aLimits.aFreqMin,aLimits.aFreqMax,aLimits.aVelMin,aLimits.aVelMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CG','F'),aT2);
    if(aMode == 1 | aMode == 3)
        hold on;    plot(aFreq1,aModel1(:,4),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
%        legend('Estimated velocity','Theor. velocity');
    end;
    if(aMode == 2 | aMode == 3)
        hold on;    plot(aFreq2,aModel2(:,4),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
%        legend('Estimated velocity','Theor. velocity');
    end;
    AtnX1 = 0.7;
    AtnX2 = 0.25;
end;  
gwlPlotFunction(aFreqDef,aModelDef(:,5),AtnX1,ASy,AtnX2,0.22,aLimits.aFreqMin,aLimits.aFreqMax,0,aLimits.aAtnMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),aT3);
if(aMode == 1 | aMode == 3)
    hold on;    plot(aFreq1,aModel1(:,5),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
%    legend('Estimated atn.','Theor. atn.');
end;
if(aMode == 2 | aMode == 3)
    hold on;    plot(aFreq2,aModel2(:,5),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
%    legend('Estimated atn.','Theor. atn.');
end;
