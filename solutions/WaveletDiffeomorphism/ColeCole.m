function ColeCole()
% In this example consider an isotropic linear viscoelastic mudium, which can be 
% described using Cole-Cole model (Jian-Fei Lu and Andrzej Hanyga, Numerical 
% modelling method for wave propagation in a linear viscoelastic medium with 
% singular memory, Geophys. J. Int. 159(2004), p. 688-702). For this mudium we 
% demonstrate the wavelet propagator for signal and its autocorrelation.
% 
% [1] M.A.Kulesh, M.Holschneider, M.S.Diallo, K.Kurennaya and F.Scherbaum Modeling 
%     of wave dispersion using continuous wavelet transforms: incorporating causality 
%     constraint with non-linear frequency-dependent attenuation // Eos Trans. AGU, 
%     86(52), Fall Meet. Suppl., Abstract S33A-0289 (2005).
% 
% FIGURE 1. (a) Phase and group velocities and (b) attenuation curve for Cole-
% Cole model.
% 
% FIGURE 2. Synthetic example: (a) propagated Shanon wavelet using Fourier method 
% (red) and the wavelet operator (black dashed line), (b) modulus of the wavelet 
% transforms of the different waveforms and (c) the corresponding phase pictures.
% 
% FIGURE 3. Synthetic example: (a) propagated cross-correlation using Fourier 
% method (red) and the wavelet operator (black dashed line), (b) modulus of the 
% wavelet transforms of the different waveforms and (c) the corresponding phase 
% pictures.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aModelName = 'model.dat';
aSignalName = 'signal.dat';
aSignalPropName = 'signalprop.dat';
aSpectrName = 'spectrum.dat';

%---------------------------------------------------------------------------
aFreq = gwlCreateAxis(256,0.0001,35,'lin',aFreqName,'Frequency');
[aFreq, aModel] = gwlDispModel(aFreqName, 'colecole', '7.87E+06,0.4,4.73E-04,1.717E-04', 'colecole', '0', aModelName);

aTime = gwlCreateAxis(1024,0,5.1175,'lin',aTimeName,'Time');
gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=shanon --wavpar=1.3 --time=0.5 --freq=8 --outtype=1 --outfile=' aSignalName])
[aTime,aSignal] = gwlSignalRead(1,aSignalName,'func',['--format=ASCII --istime --mult=0.39 --nomess'],aSignalName,'Shanon wavelet');

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aFreq,aModel(:,3),0.07,0.3,0.4,0.4,min(aFreq),max(aFreq),2800,3000,gwlGetNotation('FREQ'),gwlGetNotation('DISP','VEL'),'(a)');
    hold on;    plot(aFreq,aModel(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend(gwlGetNotation('DISP','CP','F'),gwlGetNotation('DISP','CG','F'));
gwlPlotFunction(aFreq,aModel(:,5),0.57,0.3,0.4,0.4,min(aFreq),max(aFreq),0,0.0015,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),'(b)');

%---------------------------------------------------------------------------
figure(2);
gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalPropName ' --model=' aModelName ' --step=5 --dist=2300']);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);
gwlCwt(1, aSignalName, aFreqName, 1, 'cauchy', 10, aSpectrName,'wavelet spectrum before diffeomorphism');
gwlExec('gwlDiffeoDisp',[' --infile=' aSpectrName ' --outfile=' aSpectrName ' --model=' aModelName ' --prop=3 --step=5 --dist=2300']);
[aTimeProp2,aSignalProp2] = gwlIwt(1, aSpectrName, 'delta');
[aTime,aFreq,aCwtAbs] = gwlConvert('3','',aSpectrName);
[aTime,aFreq,aCwtArg] = gwlConvert('5','--filter=5',aSpectrName);
aCwtAbsSumm = aCwtAbs(:,:,1)+aCwtAbs(:,:,2)+aCwtAbs(:,:,3)+aCwtAbs(:,:,4)+aCwtAbs(:,:,5)+aCwtAbs(:,:,6);
aCwtArgSumm = aCwtArg(:,:,1)+aCwtArg(:,:,2)+aCwtArg(:,:,3)+aCwtArg(:,:,4)+aCwtArg(:,:,5)+aCwtArg(:,:,6);
gwlPlotSeis(aTimeProp1,aSignalProp1,0.07,0.69,0.9,0.3,min(aTimeProp1),max(aTimeProp1),1,'',gwlGetNotation('TRN'),0,'(a)');
    gwlPlotSeisAdd(aTimeProp2,aSignalProp2,1,0,1);
    grid off;
gwlPlotImage(aTime,aFreq,aCwtAbsSumm,0.07,0.38,0.9,0.3,'',gwlGetNotation('FREQ'),'(b)');
gwlPlotImage(aTime,aFreq,aCwtArgSumm,0.07,0.07,0.9,0.3,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(c)');

%---------------------------------------------------------------------------
figure(3);
gwlExec('gwlAutoCorr',[' --infile=' aSignalName ' --outfile=' aSignalName]);
gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalPropName ' --model=' aModelName ' --step=5 --dist=2300 --acorr']);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);
gwlCwt(1, aSignalName, aFreqName, 1, 'cauchy', 10, aSpectrName,'wavelet spectrum before diffeomorphism');
gwlExec('gwlDiffeoDisp',[' --infile=' aSpectrName ' --outfile=' aSpectrName ' --model=' aModelName ' --prop=3 --step=5 --dist=2300 --acorr']);
[aTimeProp2,aSignalProp2] = gwlIwt(1, aSpectrName, 'delta');
[aTime,aFreq,aCwtAbs] = gwlConvert('3','',aSpectrName);
[aTime,aFreq,aCwtArg] = gwlConvert('5','--filter=5',aSpectrName);
aCwtAbsSumm = aCwtAbs(:,:,1)+aCwtAbs(:,:,2)+aCwtAbs(:,:,3)+aCwtAbs(:,:,4)+aCwtAbs(:,:,5)+aCwtAbs(:,:,6);
aCwtArgSumm = aCwtArg(:,:,1)+aCwtArg(:,:,2)+aCwtArg(:,:,3)+aCwtArg(:,:,4)+aCwtArg(:,:,5)+aCwtArg(:,:,6);
gwlPlotSeis(aTimeProp1,aSignalProp1,0.07,0.69,0.9,0.3,min(aTimeProp1),max(aTimeProp1),50,'',gwlGetNotation('TRN'),0,'(a)');
    gwlPlotSeisAdd(aTimeProp2,aSignalProp2,50,0,1);
    grid off;
gwlPlotImage(aTime,aFreq,aCwtAbsSumm,0.07,0.38,0.9,0.3,'',gwlGetNotation('FREQ'),'(b)');
gwlPlotImage(aTime,aFreq,aCwtArgSumm,0.07,0.07,0.9,0.3,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(c)');

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName);  delete(aFreqName);  delete(aModelName);  delete(aSignalName); delete(aSignalPropName);  delete(aSpectrName); 
clear all;

print -f1 -r600 -depsc ColeColeFig1;
print -f2 -r600 -depsc ColeColeFig2;
print -f3 -r600 -depsc ColeColeFig3;
