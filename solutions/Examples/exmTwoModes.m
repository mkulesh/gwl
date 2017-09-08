function exmTwoModes()

path(path, '../../mshell');
aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSignalPropName = 'signalprop.dat';
aFmax = 300;
aYmax = 50;
aDist = 113;

gwlCreateAxis(2048,0,1,'lin',aTimeName,'Time');
gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=cauchy --wavpar=5 --time=0.2 --freq=70 --outtype=1 --outfile=' aSignalName])
[aTime,aSignal,aPar] = gwlSignalRead(1,aSignalName,'func',['--format=ASCII --istime --nomess'],aSignalName,'Morlet wavelet');
[aFreqFT,aFour] = gwlFourTrans('mat',aSignal,2,aPar.aSample);

aFreq = gwlCreateAxis(1024,9.7704E-4,2047,'lin',aFreqName,'Frequency');
[aFreq, aMode1] = gwlDispModel(aFreqName, 'vel --analyt', '1256,-62,50', 'polin', '0', 'mode1.dat');
[aFreq, aMode2] = gwlDispModel(aFreqName, 'vel --analyt', '1275,125,150', 'polin', '0', 'mode2.dat');

gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=sigmode1.dat --model=mode1.dat --step=6 --decrs=2 --dist=' num2str(aDist) ]);
gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=sigmode2.dat --model=mode2.dat --step=6 --decrs=2 --dist=' num2str(aDist) ]);
gwlExec('gwlSignalSum',[' --infile=sigmode1.dat,sigmode2.dat --outfile=' aSignalPropName ]);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aFreq,aMode1(:,3),0.07,0.6,0.25,0.25,0,aFmax,1190,1445,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CP','F'),'(a)');
    hold on;    plot(aFreq,aMode2(:,3),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend('Mode 1','Mode 2');
gwlPlotFunction(aFreq,aMode1(:,4),0.40,0.6,0.25,0.25,0,aFmax,1190,1445,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CG','F'),'(b)');
    hold on;    plot(aFreq,aMode2(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend('Mode 1','Mode 2');
gwlPlotFunction(aFreqFT,abs(aFour),0.72,0.6,0.25,0.25,0,aFmax,0,max(abs(aFour)),gwlGetNotation('FREQ'),'','(c)');
gwlPlotSeis(aTimeProp1,aSignalProp1,0.07,0.13,0.9,0.4,min(aTimeProp1),max(aTimeProp1),aYmax,gwlGetNotation('TIME'),gwlGetNotation('TRN'),0,'(d)');
    grid off;

%---------------------------------------------------------------------------
figure(2);
[aFreq, aMode3] = gwlDispModel(aFreqName, 'twogauss --analyt', ['1256,-62,50,1275,125,150,' num2str(aDist)], 'polin', '0', 'mode3.dat');
gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalPropName ' --model=mode3.dat --step=6 --dist=' num2str(aDist) ]);
fid = fopen(aSignalPropName,'r'); [aTimeProp2,aSignalProp2]=gwlReadSignal(fid); fclose(fid);

gwlPlotFunction(aFreq,aMode3(:,1),0.07,0.6,0.25,0.25,0,aFmax,0,0.25,gwlGetNotation('FREQ'),gwlGetNotation('DISP','WN','F'),'(a)');
gwlPlotFunction(aFreq,aMode3(:,3),0.40,0.6,0.25,0.25,0,aFmax,1190,1445,gwlGetNotation('FREQ'),gwlGetNotation('DISP','VEL'),'(b)');
    hold on;    plot(aFreq,aMode3(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend('Phase velocity','Group velocity');
gwlPlotFunction(aFreq,aMode3(:,5),0.72,0.6,0.25,0.25,0,aFmax,-0.02,0.06,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),'(c)');
gwlPlotSeis(aTimeProp1,aSignalProp1,0.07,0.13,0.9,0.4,min(aTimeProp1),max(aTimeProp1),aYmax,gwlGetNotation('TIME'),gwlGetNotation('TRN'),0,'(d)');
    gwlPlotSeisAdd(aTimeProp2,aSignalProp2,aYmax,0,1);    
    grid off;
    
%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete('mode1.dat'); delete('mode2.dat'); delete('mode3.dat'); 
delete('sigmode1.dat'); delete('sigmode2.dat'); delete(aSignalName); delete(aSignalPropName);
clear all;
    
print -f1 -r600 -depsc exmTwoModes;
