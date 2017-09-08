function exmDiffeoLin()
% FIGURE 1.
% We choose a propagated dispersive Ricker wavelet as signal and a Morlet 
% wavelet for analysis and reconstruction. Then the linear diffeomorphism
% is compensated by dividing by a complex constant and we can simply 
% reconstruct the given signal which is the same as source signal. 
% [1] Q.Xie, M.Holschneider, M.Kulesh
%     Some remarks on linear diffeomorphisms in wavelet space // Preprint 
%     Series DFG SPP 1114, University of Bremen. Preprint 37 (2003).
%
% FIGURE 2.
% Linear diffeomorphism applied to a complex Morlet wavelet. 
% The frequency content and position of the signal is unaltered as 
% shown from the signals and spectra comparsion.
% [2] M.S.Diallo, M.Holschneider, M.Kulesh, F.Scherbaum and F.Adler
%     Characterization of the Rayleigh wave polarization attributes with 
%     continuous wavelet transform // Geophysical Research Abstracts, 
%     Vol. 5, 11237 (2003).

%---------------------------------------------------------------------------
path(path, '../../mshell');

aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';
aDiffeoName = 'cwtdiff.dat';

%---------------------------------------------------------------------------
figure(1);

gwlCreateAxis(1024,0,0.999023,'lin',aTimeName,'Time');
gwlCreateAxis(256,1,100,'lin',aFreqName,'Frequency');

aYmax = 1;
[aTime,aSignal,aPar] = gwlSignalGen(2,aTimeName,'rickdiss','200,0,800,1300,500,30',aSignalName,'Propagated Ricker wavelet');

gwlCwt(2, aSignalName, aFreqName, 1, 'morlet', 1, aSpectrName,'wavelet spectrum before diffeomorphism');

[aTime,aFreq,aCwt,aParams] = gwlConvert('3,5','--filter=1',aSpectrName);
gwlPlotImage(aTime, aFreq, aCwt(:,:,1), 0.07,0.35,0.45,0.25, '', gwlGetNotation('FREQ'), '(c)',aParams.aName);
gwlPlotImage(aTime, aFreq, aCwt(:,:,2), 0.07,0.08,0.45,0.25, gwlGetNotation('TIME'), gwlGetNotation('FREQ'), '(d)');

gwlExec('gwlDiffeoLin',[' --infile=' aSpectrName ' --outfile=' aDiffeoName ' --diffpar=1,0,3,1,0,0 --wavelet=morlet --wavpar=1']);

[aTime,aFreq,aCwt,aParams] = gwlConvert('3,5','--filter=1',aDiffeoName);
gwlPlotImage(aTime, aFreq, aCwt(:,:,1), 0.53,0.35,0.45,0.25, '', '', '(e)',aParams.aName);
gwlPlotImage(aTime, aFreq, aCwt(:,:,2), 0.53,0.08,0.45,0.25, gwlGetNotation('TIME'), '', '(f)');

[aTimeInv,aSignalInv,aParInv] = gwlIwt(2, aDiffeoName, 'morlet', 1, 1);
gwlPlotFunction(aTime,real(aSignal),0.07,0.70,0.6,0.28,min(aTime),max(aTime),-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('SIG','T'),'(a)');
    hold on; plot(aTimeInv,real(aSignalInv),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;
    legend('before diffeomorphism', 'after diffeomorphism');

[aFreqFour,aFour] = gwlFourTrans('mat',aSignal,2,aPar.aSample);
[aFreqInv,aFourInv] = gwlFourTrans('mat',aSignalInv,2,aParInv.aSample);
gwlPlotFunction(aFreqFour,abs(aFour),0.7,0.7,0.28,0.28,min(aFreq),max(aFreq),0,max(abs(aFour)),gwlGetNotation('FREQ'),'','(b)');
    hold on; plot(aFreqInv,abs(aFourInv),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;
    legend('before diffeomorphism', 'after diffeomorphism');
    
%---------------------------------------------------------------------------
figure(2);

gwlCreateAxis(1024,0,5.115,'lin',aTimeName,'Time');
gwlCreateAxis(256,1,15,'lin',aFreqName,'Frequency');
aYmax = 10;
aFreq0 = 7;

gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=morlet --wavpar=1 --time=2.5 --freq=' num2str(aFreq0) ' --outtype=1 --outfile=' aSignalName])
[aTime,aSignal,aPar] = gwlSignalRead(2,aSignalName,'func',['--format=ASCII --istime --nomess'],aSignalName,'Morlet wavelet');

gwlCwt(2, aSignalName, aFreqName, 1, 'morlet', 1, aSpectrName,'wavelet spectrum before diffeomorphism');
[aTime,aFreq,aCwt,aParams] = gwlConvert('3,5','--filter=1',aSpectrName);
gwlPlotFunction(aTime,real(aSignal),0.07,0.78,0.45,0.1,min(aTime),max(aTime),-aYmax,aYmax,'',' ',gwlGetNotation('CSIG','T',1));
    hold on; plot(aTime,imag(aSignal),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;
    set(gca,'Visible','Off');
gwlPlotImage(aTime, aFreq, aCwt(:,:,1), 0.07,0.475,0.45,0.3, '', gwlGetNotation('FREQ'), '(a)');
   gwlText(3.1,10,gwlGetNotation('CSIG','WABS',1));
gwlPlotImage(aTime, aFreq, aCwt(:,:,2), 0.07,0.17,0.45,0.3, gwlGetNotation('TIME'), gwlGetNotation('FREQ'), '(b)');
   gwlText(3.1,10,gwlGetNotation('CSIG','WARG',1));
    
gwlExec('gwlDiffeoLin',[' --infile=' aSpectrName ' --outfile=' aDiffeoName ' --diffpar=1,0,3.3,1,0,0 --wavelet=cauchy --wavpar=6']);
[aTime,aFreq,aCwt,aParams] = gwlConvert('3,5','--filter=1',aDiffeoName);
[aTimeInv,aSignalInv,aParInv] = gwlIwt(2, aDiffeoName, 'cauchy', 5, 1);
gwlPlotFunction(aTime,real(aSignalInv),0.54,0.78,0.45,0.1,min(aTime),max(aTime),-aYmax,aYmax,'','',gwlGetNotation('CSIG','T',2));
    hold on; plot(aTime,imag(aSignalInv),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;
    set(gca,'Visible','Off');
gwlPlotImage(aTime, aFreq, aCwt(:,:,1), 0.54,0.475,0.45,0.3, '', ' ', '(c)');
   gwlText(2.7,10,gwlGetNotation('CSIG','WABS',2));
gwlPlotImage(aTime, aFreq, aCwt(:,:,2), 0.54,0.17,0.45,0.3, gwlGetNotation('TIME'), ' ', '(d)');
   gwlText(2.7,10,gwlGetNotation('CSIG','WARG',2));

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aDiffeoName); 
clear all;

print -f1 -r600 -depsc exmDiffeoLinFig1;
print -f2 -r600 -depsc exmDiffeoLinFig2;

