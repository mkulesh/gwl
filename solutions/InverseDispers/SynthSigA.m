function SynthSigA()
% In this example we compare of the dispersion models formulated in the Fourier 
% and wavelet spaces for complex synthetic signal. As an approximation of the 
% phase velocity function, we use the three-parameter exponential approximation.
% 
% [1] M.A.Kulesh, M.S.Diallo and M.Holschneider Wavelet analysis of ellipticity, 
%     dispersion, and dissipation properties of Rayleigh waves // Acoustical Physics. 
%     V. 51. No. 4. P. 421-434 (2005).

%---------------------------------------------------------------------------
path(path, '../../mshell');
aFreqName = 'freq.dat';
aModelName = 'model.dat';
aSignalName = 'signal.dat';
aSignalPropName = 'signalprop.dat';
aSpectrName = 'spectrum.dat';
aSpectrOptName = 'spectrumopt.dat';
aModelOpt1Name = 'modelopt1.dat';
aModelOpt2Name = 'modelopt2.dat';
aTmin = 0.0;
aTmax = 5;
aYmax = 0.05;
aInd = 1001;

%---------------------------------------------------------------------------
aFreq = gwlCreateAxis(128,0.1,20,'lin',aFreqName,'Frequency');
[aFreq, aModel] = gwlDispModel(aFreqName, 'vel', '1300,300,10 --analyt', 'polin', '0',aModelName);
[aTime, aSignal] = gwlSignalRead(2,'SynthSigA.asc','func','--istime',aSignalName,'Synthetic complex signal');

gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalPropName ' --model=' aModelName ' --step=1 --dist=2000']);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);

gwlDispModel(aFreqName, 'gauss', '0.0007,-0.00001,5 --analyt --nomess', 'polin', '0', aModelName, 'Initial dispersion model');
aFreq = gwlCreateAxis(128,0.1,20,'lin --sign=full',aFreqName,'Frequency');
gwlCwt(2, aSignalPropName, aFreqName, 2, 'morlet', 1, aSpectrName,'Wavelet spectrum');
gwlExec('gwlOptiSP',[' --infile=' aModelName ' --outfile=' aModelOpt1Name ' --spec=' aSpectrName ' --dist=2000 --cmpl=3 --prog --name="modulus optimized model"']);
gwlExec('gwlOptiSP',[' --infile=' aModelOpt1Name ' --outfile=' aModelOpt2Name ' --spec=' aSpectrName ' --ospec=' aSpectrOptName ' --dist=2000 --cmpl=4 --prog --name="argument optimized model"']);
[aTimeProp2,aSignalProp2] = gwlIwt(2, aSpectrOptName, 'delta');

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aTimeProp1,real(aSignalProp1(:,1)),0.07,0.84,0.91,0.13,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('CSIG','T',1),'(a)');
    hold on;    plot(aTimeProp1,imag(aSignalProp1(:,1)),'--black','LineWidth',1);    hold off;
    
gwlPlotFunction(aTimeProp1,real(aSignalProp1(:,2)),0.07,0.71,0.91,0.13,aTmin,aTmax,-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('CSIG','T',2),'(b)');
    hold on;    plot(aTimeProp1,imag(aSignalProp1(:,2)),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    hold on;    plot(aTimeProp2,real(aSignalProp2(:,2)),'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
    hold on;    plot(aTimeProp2,imag(aSignalProp2(:,2)),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;

[aTimeWT,aFreqWT,aCwtAbs] = gwlConvert('3','',aSpectrName);
[aTimeWT,aFreqWT,aCwtArg] = gwlConvert('5','--filter=5',aSpectrName);
[aTimeWT,aFreqWT,aCwtOptAbs] = gwlConvert('3','',aSpectrOptName);
[aTimeWT,aFreqWT,aCwtOptArg] = gwlConvert('5','--filter=5',aSpectrOptName);

gwlPlotImage(aTimeWT(1:aInd),aFreqWT,aCwtAbs(:,(1:aInd),1),0.07,0.36,0.3,0.29,'',gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('CSIG','WABS',1)]);
    line([0,aTmax],[0,0],'Color','black');

gwlPlotImage(aTimeWT(1:aInd),aFreqWT,aCwtAbs(:,(1:aInd),2),0.375,0.36,0.3,0.29,'','',gwlGetNotation('CSIG','WABS',2));
    line([0,aTmax],[0,0],'Color','black');

gwlPlotImage(aTimeWT(1:aInd),aFreqWT,aCwtOptAbs(:,(1:aInd),2),0.68,0.36,0.3,0.29,'','',gwlGetNotation('CSIG','WABS','n'));
    line([0,aTmax],[0,0],'Color','black');

gwlPlotImage(aTimeWT(1:aInd),aFreqWT,aCwtArg(:,(1:aInd),1),0.07,0.07,0.3,0.28,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(d) ' gwlGetNotation('CSIG','WARG',1)]);
    line([0,aTmax],[0,0],'Color','black');

gwlPlotImage(aTimeWT(1:aInd),aFreqWT,aCwtArg(:,(1:aInd),2),0.375,0.07,0.3,0.28,gwlGetNotation('TIME'),'',gwlGetNotation('CSIG','WARG',2));
    line([0,aTmax],[0,0],'Color','black');

gwlPlotImage(aTimeWT(1:aInd),aFreqWT,aCwtOptArg(:,(1:aInd),2),0.68,0.07,0.3,0.28,gwlGetNotation('TIME'),'',gwlGetNotation('CSIG','WARG','n'));
    line([0,aTmax],[0,0],'Color','black');

%---------------------------------------------------------------------------
figure(2);
fid = fopen(aModelOpt1Name,'r'); [aFreq, aModelOpt1]=gwlReadDispModel(fid); fclose(fid);
fid = fopen(aModelOpt2Name,'r'); [aFreq, aModelOpt2]=gwlReadDispModel(fid); fclose(fid);

gwlPlotFunction(aFreq, aModel(:,3),0.07,0.3,0.4,0.3,0,max(aFreq),1200,1700,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CP','F'),'(e)');
    hold on;    plot(aFreq, aModelOpt1(:,3),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    hold on;    plot(aFreq, aModelOpt2(:,3),'Color',gwlGetColor(0),'LineStyle','-.','LineWidth',1);    hold off;

gwlPlotFunction(aFreq, aModel(:,4),0.55,0.3,0.4,0.3,0,max(aFreq),1100,1700,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CG','F'),'(f)');
    hold on;    plot(aFreq, aModelOpt1(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    hold on;    plot(aFreq, aModelOpt2(:,4),'Color',gwlGetColor(0),'LineStyle','-.','LineWidth',1);    hold off;

%---------------------------------------------------------------------------
pause(0.00001);
delete(aFreqName);  delete(aModelName);  delete(aSignalName);  delete(aSignalPropName);  delete(aSpectrName); 
delete(aSpectrOptName); delete(aModelOpt1Name); delete(aModelOpt2Name);
clear all;

print -f1 -r600 -depsc SynthSigAFig1;
print -f2 -r600 -depsc SynthSigAFig2;
