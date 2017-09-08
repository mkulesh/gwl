function exmColeCole()

%---------------------------------------------------------------------------
path(path, '../../mshell');
aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aModelName = 'model.dat';
aSignalName = 'signal.dat';

%---------------------------------------------------------------------------
aFreq = gwlCreateAxis(256,0.0001,35,'lin',aFreqName,'Frequency');
[aFreq, aModel] = gwlDispModel(aFreqName, 'colecole', '7.87E+06,0.4,4.73E-04,1.717E-04', 'colecole', '0', aModelName);

aTime = gwlCreateAxis(1024,0,5.1175,'lin',aTimeName,'Time');
gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=shanon --wavpar=1.3 --time=0.5 --freq=8 --outtype=1 --outfile=' aSignalName])
[aTime,aSignal] = gwlSignalRead(1,aSignalName,'func',['--format=ASCII --istime --mult=0.39 --nomess'],aSignalName,'Shanon wavelet');

gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalName ' --model=' aModelName ' --step=5 --dist=2300']);
fid = fopen(aSignalName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aFreq,aModel(:,3),0.07,0.7,0.4,0.25,min(aFreq),max(aFreq),2800,3000,gwlGetNotation('FREQ'),gwlGetNotation('DISP','VEL'),'(a)');
    hold on;    plot(aFreq,aModel(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend(gwlGetNotation('DISP','CP','F'),gwlGetNotation('DISP','CG','F'));
gwlPlotFunction(aFreq,aModel(:,5),0.57,0.7,0.4,0.25,min(aFreq),max(aFreq),0,0.0015,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),'(b)');
gwlPlotSeis(aTimeProp1,aSignalProp1,0.07,0.07,0.9,0.55,min(aTimeProp1),max(aTimeProp1),1,gwlGetNotation('TIME'),gwlGetNotation('SIG','T','m'),0,'(c)');

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName);  delete(aFreqName);  delete(aModelName);  delete(aSignalName);
clear all;

print -f1 -r600 -depsc exmColeCole;
