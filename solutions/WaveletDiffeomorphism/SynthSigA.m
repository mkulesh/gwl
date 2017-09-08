function SynthSigA()
% To test wavelet-based propagation operator which can properly handle dispersion 
% and dissipation, we simulate the propagation of a Ricker (RICKER, 1935) waveform 
% (central frequency around 30 Hz) in a dispersive and attenuating medium. The 
% signals propagated at four successive stations using the Fourier integral method 
% and wavelet diffeomorphism are shown in Figure 1d for comparison. Note the 
% excellent agreement between the results from Fourier and wavelet propagation. 
% Figures 2a,b show the wavelet modulus and the phase images for the propagated 
% waveforms. Considering the deformation pattern of modulus and phase images, we 
% see that the propagation operator manipulates the phase and the modulus in 
% different ways, as elucidated in the mathematical treatment in [2].
% 
% [1] M.Kulesh, M.Holschneider, M.S.Diallo, Q.Xie and F.Scherbaum Modeling of wave 
%     dispersion using continuous wavelet transforms // Preprint Series DFG SPP 1114, 
%     University of Bremen. Preprint 40 (2003).
% 
% [2] M.Kulesh, M.Holschneider, M.S.Diallo, Q.Xie and F.Scherbaum Modeling of wave 
%     dispersion using continuous wavelet transforms // Pure and Applied Geophysics. 
%     V. 162. No. 5. P. 843-855 (2005).
% 
% [3] M.S.Diallo, M.Holschneider, F.Scherbaum and M.Kulesh Characterization of 
%     dispersive Rayleigh waves using wavelet transform // Eos. Trans. AGU, 84(46), 
%     Fall Meet. Suppl., Abstract S22B-0442 (2003).
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
% FIGURE 1. (a) Wavenumber curves, (b) phase and group velocity, (c) attenuation 
% coefficient, (d) propagated Ricker wavelet using Fourier method (red) and the 
% wavelet operator (black line).
% 
% FIGURE 2. (a) Modulus of the wavelet transforms of the different waveforms in 
% (Figure 1d) obtained with wavelet method and (b) the corresponding phase 
% pictures. Note the difference in the deformation of the phase and the modulus of 
% the transforms.
% 
% FIGURE 3. The same as Fig. 1-2, was used in [4].
% 
% FIGURE 4. Synthetic example: (a) Comparison between original cross-correlations 
% Tij(t) (red) and those obtained with the wavelet propagator (blue), (b) wavelet-
% modulus of the different cross-correlations in (a), and (c) the deformed 
% wavelet-phase pictures of the crosscorrelations in (a).


%---------------------------------------------------------------------------
path(path, '../../mshell');
aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aModelName = 'model.dat';
aSignalName = 'signal.dat';
aCorrName = 'autocorr.dat';
aSignalPropName = 'signalprop.dat';
aSpectrName = 'spectrum.dat';
aSpectrPropName = 'spectrumprop.dat';
aYmax = 1.5;

%---------------------------------------------------------------------------
aTime = gwlCreateAxis(1024,0,2.046,'lin',aTimeName,'Time');
aFreq = gwlCreateAxis(256,1,100,'lin',aFreqName,'Frequency');
[aFreq, aModel] = gwlDispModel(aFreqName, 'vel', '1300,200,30', 'polin', '+3.602E-04,+1.849E-05,-1.600E-06,+3.000E-08,-1.668E-10',aModelName);
[aTime, aSignal] = gwlSignalGen(1,aTimeName,'rickdiss','200,0,500,1300,200,30',aSignalName,'Propagated Ricker wavelet');

%---------------------------------------------------------------------------
figure(1);
gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalPropName ' --model=' aModelName ' --step=3 --dist=500']);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);
gwlCwt(1, aSignalName, aFreqName, 1, 'morlet', 2, aSpectrName,'wavelet spectrum before diffeomorphism');
gwlExec('gwlDiffeoDisp',[' --infile=' aSpectrName ' --outfile=' aSpectrPropName ' --model=' aModelName ' --step=3 --dist=500']);
[aTimeProp2,aSignalProp2] = gwlIwt(1, aSpectrPropName, 'delta');
[aTime,aFreq,aCwtAbs] = gwlConvert('3','',aSpectrPropName);
[aTime,aFreq,aCwtArg] = gwlConvert('5','--filter=20',aSpectrPropName);
gwlPlotFunction(aFreq,aModel(:,1),0.07,0.7,0.26,0.25,min(aFreq),max(aFreq),0,0.1,gwlGetNotation('FREQ'),gwlGetNotation('DISP','WN','F'),'(a)');
gwlPlotFunction(aFreq,aModel(:,3),0.39,0.7,0.26,0.25,min(aFreq),max(aFreq),1200,1600,gwlGetNotation('FREQ'),gwlGetNotation('DISP','VEL'),'(b)');
    hold on;    plot(aFreq,aModel(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend(gwlGetNotation('DISP','CP','F'),gwlGetNotation('DISP','CG','F'));
gwlPlotFunction(aFreq,aModel(:,5),0.71,0.7,0.26,0.25,min(aFreq),max(aFreq),0,0.0006,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),'(c)');
gwlPlotSeis(aTimeProp1,aSignalProp1,0.07,0.1,0.9,0.5,min(aTimeProp1),max(aTimeProp1),aYmax,gwlGetNotation('TIME'),gwlGetNotation('TRN'),0,'(d)');
    gwlPlotSeisAdd(aTimeProp2,aSignalProp2,aYmax,0,1);
    grid off;

%---------------------------------------------------------------------------
figure(2);
gwlPlotImage(aTime,aFreq,aCwtAbs(:,:,1)+aCwtAbs(:,:,2)+aCwtAbs(:,:,3)+aCwtAbs(:,:,4),0.07,0.53,0.9,0.42,'',gwlGetNotation('FREQ'),'(a)');
    gwlText(0.3,85, gwlGetNotation('SIG','WABS','1'));
    gwlText(0.67,85,gwlGetNotation('SIG','WABS','2'));
    gwlText(1.05,85,gwlGetNotation('SIG','WABS','3'));
    gwlText(1.4,85, gwlGetNotation('SIG','WABS','4'));
gwlPlotImage(aTime,aFreq,aCwtArg(:,:,1)+aCwtArg(:,:,2)+aCwtArg(:,:,3)+aCwtArg(:,:,4),0.07,0.1,0.9,0.41,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(b)');
    gwlText(0.3,85, gwlGetNotation('SIG','WARG','1'));
    gwlText(0.67,85,gwlGetNotation('SIG','WARG','2'));
    gwlText(1.05,85,gwlGetNotation('SIG','WARG','3'));
    gwlText(1.4,85, gwlGetNotation('SIG','WARG','4'));

%---------------------------------------------------------------------------
figure(3);
gwlPlotFunction(aFreq,aModel(:,1),0.07,0.77,0.26,0.2,min(aFreq),max(aFreq),0,0.1,gwlGetNotation('FREQ'),gwlGetNotation('DISP','WN','F'),'(a)');
gwlPlotFunction(aFreq,aModel(:,3),0.39,0.77,0.26,0.2,min(aFreq),max(aFreq),1200,1600,gwlGetNotation('FREQ'),gwlGetNotation('DISP','VEL'),'(b)');
    hold on;    plot(aFreq,aModel(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend(gwlGetNotation('DISP','CP','F'),gwlGetNotation('DISP','CG','F'));
gwlPlotFunction(aFreq,aModel(:,5),0.71,0.77,0.26,0.2,min(aFreq),max(aFreq),0,0.0006,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),'(c)');
aSig1(:,1) = aSignalProp1(:,1);
aSig1(:,2) = aSignalProp1(:,4);
aSig2(:,1) = aSignalProp2(:,1);
aSig2(:,2) = aSignalProp2(:,4);
gwlPlotSeis(aTimeProp1,aSig1,0.07,0.5,0.9,0.2,min(aTimeProp1),max(aTimeProp1),aYmax,gwlGetNotation('TIME'),gwlGetNotation('TRN'),0,'(d)');
    gwlPlotSeisAdd(aTimeProp2,aSig2,aYmax,0,1);
    grid off;
gwlPlotImage(aTime,aFreq,aCwtAbs(:,:,1)+aCwtAbs(:,:,4),0.07,0.29,0.9,0.2,'',gwlGetNotation('FREQ'),'(e)');
    gwlText(0.4,85,gwlGetNotation('SIG','WABS','1'));
    gwlText(1.5,85,gwlGetNotation('SIG','WABS','2'));
gwlPlotImage(aTime,aFreq,aCwtArg(:,:,1)+aCwtArg(:,:,4),0.07,0.08,0.9,0.2,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(f)');
    gwlText(0.4,85,gwlGetNotation('SIG','WARG','1'));
    gwlText(1.5,85,gwlGetNotation('SIG','WARG','2'));

%---------------------------------------------------------------------------
figure(4);
gwlExec('gwlAutoCorr',[' --infile=' aSignalName ' --outfile=' aCorrName]);
gwlExec('gwlDiffeoDisp',[' --infile=' aCorrName ' --outfile=' aSignalPropName ' --model=' aModelName ' --step=3 --dist=500 --acorr']);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);
gwlCwt(1, aCorrName, aFreqName, 1, 'morlet', 2, aSpectrName,'wavelet spectrum before diffeomorphism');
gwlExec('gwlDiffeoDisp',[' --infile=' aSpectrName ' --outfile=' aSpectrPropName ' --model=' aModelName ' --step=3 --dist=500 --acorr']);
[aTimeProp2,aSignalProp2] = gwlIwt(1, aSpectrPropName, 'delta');
[aTime,aFreq,aCwtAbs] = gwlConvert('3','',aSpectrPropName);
[aTime,aFreq,aCwtArg] = gwlConvert('5','--filter=20',aSpectrPropName);
gwlPlotSeis(aTimeProp1,aSignalProp1,0.07,0.68,0.9,0.3,min(aTimeProp1),max(aTimeProp1),10*aYmax,'',gwlGetNotation('TRN'),0,'(a)');
    gwlPlotSeisAdd(aTimeProp2,aSignalProp2,10*aYmax,0,1);
    gwlText(0.05,16*aYmax, gwlGetNotation('TSIG','T','1-1'));
    gwlText(0.05,26*aYmax,gwlGetNotation('TSIG','T','1-2'));
    gwlText(0.05,36*aYmax,gwlGetNotation('TSIG','T','1-3'));
    gwlText(0.05,46*aYmax,gwlGetNotation('TSIG','T','1-4'));
    grid off;
gwlPlotImage(aTime,aFreq,aCwtAbs(:,:,1)+aCwtAbs(:,:,2)+aCwtAbs(:,:,3)+aCwtAbs(:,:,4),0.07,0.37,0.9,0.3,'',gwlGetNotation('FREQ'),'(b)');
gwlPlotImage(aTime,aFreq,aCwtArg(:,:,1)+aCwtArg(:,:,2)+aCwtArg(:,:,3)+aCwtArg(:,:,4),0.07,0.06,0.9,0.3,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(c)');

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName);  delete(aFreqName);  delete(aModelName);  delete(aSignalName);
delete(aCorrName);  delete(aSignalPropName);  delete(aSpectrName);  delete(aSpectrPropName);
clear all;

print -f1 -r600 -depsc SynthSigAFig1;
print -f2 -r600 -depsc SynthSigAFig2;
print -f3 -r600 -depsc SynthSigAFig3;
print -f4 -r600 -depsc SynthSigAFig4;
