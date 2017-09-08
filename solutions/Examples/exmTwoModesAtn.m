function exmTwoModesAtn()

%---------------------------------------------------------------------------
path(path, '../../mshell');
aLimits.aVelMin = 900;
aLimits.aVelMax = 2000;
aLimits.aAtnMax = 0.002;

gwlExec('gwlCreateAxis', '--outfile=exmTwoModesAtnBTime.dat --count=1024 --min=0 --max=1.023 --name="Time" --nomess');
gwlExec('gwlCreateAxis', '--outfile=exmTwoModesAtnBFreq.dat --count=1024 --min=0.1 --max=100 --name="Frequency" --nomess');
gwlExec('gwlDispModel', '--infile=exmTwoModesAtnBFreq.dat --outfile=exmTwoModesAtnBMod1.dat --analyt --wn=vel --wnpar=1100,400,20 --atn=polin --atnpar=6.204e-04,3.698e-05,-3.200e-06,6.000e-08,-3.336e-10 --name="1th dispersion model" --nomess');
gwlExec('gwlSignalGen', '--infile=exmTwoModesAtnBTime.dat --outfile=exmTwoModesAtnBSour1.dat --type=rickdiss --par=120,0,150,1100,400,20 --name="Propagated Ricker wavelet" --nomess');
gwlExec('gwlDiffeoDisp', '--infile=exmTwoModesAtnBSour1.dat --outfile=exmTwoModesAtnBSig1.dat --model=exmTwoModesAtnBMod1.dat --step=4 --dist=170 --name="1th propagation mode"');
gwlExec('gwlDispModel', '--infile=exmTwoModesAtnBFreq.dat --outfile=exmTwoModesAtnBMod2.dat --analyt --wn=vel --wnpar=1200,800,50 --atn=polin --atnpar=0.00120,3.898e-005,-3.16e-06,6.03e-08,-3.336e-10 --name="2nd dispersion model" --nomess');
gwlExec('gwlSignalGen', '--infile=exmTwoModesAtnBTime.dat --outfile=exmTwoModesAtnBSour2.dat --type=rickdiss --par=70,60,150,1200,800,50 --name="Propagated Ricker wavelet" --nomess');
gwlExec('gwlDiffeoDisp', '--infile=exmTwoModesAtnBSour2.dat --outfile=exmTwoModesAtnBSig2.dat --model=exmTwoModesAtnBMod2.dat --step=4 --dist=170 --name="2nd propagation mode"');
gwlExec('gwlSignalSum', '--infile=exmTwoModesAtnBSig1.dat,exmTwoModesAtnBSig2.dat --outfile=exmTwoModesAtnBSig.dat --name="Synthetic seismogram with two modes" --nomess');
   
fid = fopen('exmTwoModesAtnBMod1.dat','r'); [aFreq,aModel1]=gwlReadDispModel(fid); fclose(fid);
fid = fopen('exmTwoModesAtnBMod2.dat','r'); [aFreq,aModel2]=gwlReadDispModel(fid); fclose(fid);
fid = fopen('exmTwoModesAtnBSig.dat','r'); [aTime,aSeis,aParSeis]=gwlReadSignal(fid); fclose(fid);

%---------------------------------------------------------------------------
figure(1);
aLimits.aFreqMin = 0;
aLimits.aFreqMax = 100;
aAmplMax = 0.7;
[aFreqFT,aFour] = gwlFourTrans('gwl','exmTwoModesAtnBSig.dat',2);
gwlPlotFunction(aFreq,aModel1(:,3),0.07,0.71,0.25,0.25,aLimits.aFreqMin,aLimits.aFreqMax,aLimits.aVelMin,aLimits.aVelMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CP','F'),'(a)');
    hold on;    plot(aFreq,aModel2(:,3),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend('Mode 1','Mode 2');
gwlPlotFunction(aFreq,aModel1(:,4),0.39,0.71,0.25,0.25,aLimits.aFreqMin,aLimits.aFreqMax,aLimits.aVelMin,aLimits.aVelMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CG','F'),'(b)');
    hold on;    plot(aFreq,aModel2(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend('Mode 1','Mode 2');
gwlPlotFunction(aFreq,aModel1(:,5),0.70,0.71,0.25,0.25,aLimits.aFreqMin,aLimits.aFreqMax,0,aLimits.aAtnMax,gwlGetNotation('FREQ'),gwlGetNotation('DISP','IMWN','F'),'(c)');
    hold on;    plot(aFreq,aModel2(:,5),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend('Mode 1','Mode 2');
gwlPlotSeis(aTime,aSeis,0.07,0.34,0.88,0.31,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),gwlGetNotation('SIG','T','m'),0,'(d)');
    gwlText(0.01,1*aAmplMax,gwlGetNotation('SIG','T','1'));
    gwlText(0.01,2*aAmplMax,gwlGetNotation('SIG','T','2'));
    gwlText(0.01,3*aAmplMax,gwlGetNotation('SIG','T','3'));
    gwlText(0.01,4*aAmplMax,gwlGetNotation('SIG','T','4'));
    gwlText(0.01,5*aAmplMax,gwlGetNotation('SIG','T','5'));
    grid off;
gwlPlotFunction(aFreqFT,abs(aFour),0.32,0.07,0.4,0.2,aLimits.aFreqMin,aLimits.aFreqMax,0,max(max(abs(aFour))),gwlGetNotation('FREQ'),gwlGetNotation('SIG','FM','m'),'(e)');
    
%---------------------------------------------------------------------------
figure(2);
aSig=1.6*aSeis(:,3);
aSamplfreq=1000;
[aFreq1,aFour1] = gwlFourTrans('mat',aSig,1,aSamplfreq);
[aFreq2,aFour2] = gwlFourTrans('mat',aSig,2,aSamplfreq);
gwlPlotFunction(aTime,aSig,0.07,0.55,0.9,0.25,min(aTime),max(aTime),-1,1,gwlGetNotation('TIME'),gwlGetNotation('SIG','T'),'(a)');
gwlPlotFunction(aFreq1,abs(aFour1),0.07,0.2,0.28,0.28,0.0,aSamplfreq,0,max(abs(aFour1)),gwlGetNotation('FREQ'),gwlGetNotation('SIG','FM'),'(b)');
gwlPlotFunction(aFreq2,abs(aFour2),0.38,0.2,0.28,0.28,-aSamplfreq/2,aSamplfreq/2,0,max(abs(aFour2)),gwlGetNotation('FREQ'),' ','(c)');
gwlPlotFunction(aFreq2,abs(aFour2),0.69,0.2,0.28,0.28,0,100,0,max(abs(aFour2)),gwlGetNotation('FREQ'),' ','(d)');

%---------------------------------------------------------------------------
pause(0.00001);
delete('exmTwoModesAtn*.dat'); 
print -f1 -r600 -depsc exmTwoModesAtnFig1;
print -f2 -r600 -depsc exmTwoModesAtnFig2;
