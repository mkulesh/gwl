function exmCft()
% exmCft(): An example of CFT calculation using matlab and gwl algorithms

%---------------------------------------------------------------------------
path(path, '../../mshell');

aTimeName = 'time.dat';
aSignalName = 'signal.dat';
aYmax = 10.0;

%---------------------------------------------------------------------------
figure(1);
gwlCreateAxis(1024,0,5.115,'lin',aTimeName,'Time');
[aTime,aSignal,aParSig] = gwlSignalGen(2,aTimeName,'harmrot','2.0,7.0,1.0,2.0,5.0,5.0,0.318',aSignalName,'Rotated harmonic function');

gwlPlotFunction(aTime,real(aSignal),0.07,0.7,0.9,0.25,min(aTime),max(aTime),-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('CSIG','T'),'(a)');
    hold on;    plot(aTime,imag(aSignal),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});

aShift = 2;
[aFreq,aFour] = gwlFourTrans('gwl',aSignalName,aShift);
[aFreq1,aFour1] = gwlFourTrans('mat',aSignal,aShift,aParSig.aSample);

aImax = 3000;
aFmax = 10;
gwlPlotFunction(aFreq,real(aFour),0.07,0.07,0.9,0.55,-aFmax,aFmax,-aImax,aImax,gwlGetNotation('FREQ'),gwlGetNotation('CSIG','REF'),'(b)');
    hold on;    plot(aFreq1,real(aFour1),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;


%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aSignalName);
clear all;

