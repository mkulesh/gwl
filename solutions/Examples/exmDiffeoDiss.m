function exmDiffeoDiss()

%---------------------------------------------------------------------------
path(path, '../../mshell');
aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aModelName = 'model.dat';
aSignalName = 'signal.dat';
aSignalPropName = 'signalprop.dat';
aSpectrName = 'spectrum.dat';
aElliparName = 'elli.dat';
aDist = 1100;
aTime = gwlCreateAxis(1024,0,2.046,'lin',aTimeName,'Time');

%---------------------------------------------------------------------------
figure(1);
gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=morlet --wavpar=0.5 --time=0.5 --freq=40 --outtype=1 --outfile=' aSignalName])
[aTime,aSignal,aPar] = gwlSignalRead(1,aSignalName,'func',['--format=ASCII --istime --nomess'],aSignalName,'Morlet wavelet');

aFreq = gwlCreateAxis(256,0.1,90,'lin',aFreqName,'Frequency');
[aFreq, aModel] = gwlDispModel(aFreqName, 'vel', '1300,400,30', 'polin', '0', aModelName);

gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalPropName ' --model=' aModelName ' --step=1 --dist=' num2str(aDist)]);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);

gwlCwt(1, aSignalName, aFreqName, 2, 'morlet', 3, aSpectrName,'wavelet spectrum before diffeomorphism');
gwlExec('gwlDiffeoDisp',[' --infile=' aSpectrName ' --outfile=' aSpectrName ' --model=' aModelName ' --prop=2 --step=1 --dist=' num2str(aDist)]);
[aTimeProp2,aSignalProp2] = gwlIwt(1, aSpectrName, 'delta');
[aTime,aFreq,aCwtAbs] = gwlConvert('3','',aSpectrName);
[aTime,aFreq,aCwtArg] = gwlConvert('5','--filter=2',aSpectrName);

aYmax = max(abs(aSignalProp1(:,1)));
gwlPlotFunction(aTimeProp1,aSignalProp1(:,1),0.07,0.78,0.45,0.1,min(aTime),max(aTime),-aYmax,aYmax,'','',gwlGetNotation('SIG','T',1));
    set(gca,'Visible','Off');
gwlPlotImage(aTime,aFreq,aCwtAbs(:,:,1),0.07,0.475,0.45,0.3,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('SIG','WABS',1)]);
gwlPlotImage(aTime,aFreq,aCwtArg(:,:,1),0.07,0.17,0.45,0.3,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('SIG','WARG',1)]);
gwlPlotFunction(aTimeProp1,aSignalProp1(:,2),0.54,0.78,0.45,0.1,min(aTimeProp1),max(aTimeProp1),-aYmax,aYmax,'','',gwlGetNotation('SIG','T',2));
    set(gca,'Visible','Off');
gwlPlotImage(aTime,aFreq,aCwtAbs(:,:,2),0.54,0.475,0.45,0.3,'',' ',['(c) ' gwlGetNotation('SIG','WABS',2)]);
   gwlText(1.15,85,gwlGetNotation('DISP','CG','F'));
   line(0.5+aDist./aModel(:,4),aFreq,'Color',gwlGetColor(1),'LineWidth',2);
gwlPlotImage(aTime,aFreq,aCwtArg(:,:,2),0.54,0.17,0.45,0.3,gwlGetNotation('TIME'),' ',['(d) ' gwlGetNotation('SIG','WARG',2)]);
   gwlText(1.15,85,gwlGetNotation('DISP','CP','F'));
   line(0.5+aDist./aModel(:,3),aFreq,'Color',gwlGetColor(1),'LineWidth',2);

%---------------------------------------------------------------------------
figure(2);
gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=morlet --wavpar=0.5 --time=0.5 --freq=20 --outtype=1 --outfile=' aSignalName])
[aTime,aSignal,aPar] = gwlSignalRead(2,aSignalName,'func',['--format=ASCII --istime --nomess'],aSignalName,'Morlet wavelet');

aFreq = gwlCreateAxis(256,0.1,60,'lin  --sign=full',aFreqName,'Frequency');
[aFreq, aModel] = gwlDispModel(aFreqName, 'vel', '1300,400,20', 'polin', '0', aModelName);

gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalPropName ' --model=' aModelName ' --step=1 --dist=' num2str(aDist)]);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);

gwlCwt(2, aSignalName, aFreqName, 2, 'morlet', 4, aSpectrName,'wavelet spectrum before diffeomorphism');
gwlExec('gwlDiffeoDisp',[' --infile=' aSpectrName ' --outfile=' aSpectrName ' --model=' aModelName ' --prop=5 --step=1 --dist=' num2str(aDist)]);
[aTimeProp2,aSignalProp2] = gwlIwt(2, aSpectrName, 'delta');
gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=complex --filter=1 --chan=0']);
[aTimeE,aFreqE,aCwtS1] = gwlConvert('1,14','--dergee --nomess',aElliparName);
gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=complex --filter=1 --chan=1']);
[aTimeE,aFreqE,aCwtS2] = gwlConvert('1,14','--dergee --nomess',aElliparName);

aYmax = max(abs(aSignalProp1(:,1)));
gwlPlotFunction(aTimeProp1,real(aSignalProp1(:,1)),0.07,0.78,0.45,0.1,min(aTime),max(aTime),-aYmax,aYmax,'','',gwlGetNotation('CSIG','T',1));
    hold on; plot(aTime,imag(aSignalProp1(:,1)),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;
    set(gca,'Visible','Off');
gwlPlotImage(aTimeE,aFreqE,aCwtS1(:,:,1),0.07,0.475,0.45,0.3,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('EPAR','RMAX','WG',1)]);
gwlPlotImage(aTimeE,aFreqE,aCwtS1(:,:,2),0.07,0.17,0.45,0.3,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','TILT','WG',1)]);
gwlPlotFunction(aTimeProp1,real(aSignalProp1(:,2)),0.54,0.78,0.45,0.1,min(aTimeProp1),max(aTimeProp1),-aYmax,aYmax,'','',gwlGetNotation('CSIG','T',2));
    hold on; plot(aTime,imag(aSignalProp1(:,2)),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;
    set(gca,'Visible','Off');
gwlPlotImage(aTimeE,aFreqE,aCwtS2(:,:,1),0.54,0.475,0.45,0.3,'',' ',['(c) ' gwlGetNotation('EPAR','RMAX','WG',2)]);
   gwlText(1.15,55,gwlGetNotation('DISP','CG','F'));
   line(0.5+aDist./aModel(:,4),aFreq,'Color',gwlGetColor(1),'LineWidth',2);
gwlPlotImage(aTimeE,aFreqE,aCwtS2(:,:,2),0.54,0.17,0.45,0.3,gwlGetNotation('TIME'),' ',['(d) ' gwlGetNotation('EPAR','TILT','WG',2)]);
   gwlText(1.15,55,gwlGetNotation('DISP','CP','F'));
   line(0.5+aDist./aModel(:,3),aFreq,'Color',gwlGetColor(1),'LineWidth',2);

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName);  delete(aFreqName);  delete(aModelName);  delete(aSignalName); delete(aSignalPropName);  delete(aSpectrName); delete(aElliparName);
clear all;

print -f1 -r600 -depsc exmDiffeoDiss;
print -f2 -r600 -depsc exmDiffeoDissET;
