function VibroBridge()
% VibroBridge(): Plotting of experimental data: spectral properties of 
% velocities and accelerations by a vibrational measurement at a bridge

%---------------------------------------------------------------------------
path(path, '../../mshell');
[aTime,aSignal,aPar] = gwlSignalRead(1,'VibroBridge.asc','seis','--format=ASCII --istime --tmin=0');
aSFreq = aPar.aSample;
 
%---------------------------------------------------------------------------
figure(2);
localPlotRes(aTime,aSignal(:,1),aPar.aSample,'Geophone: 1, Cannel: Z');

%---------------------------------------------------------------------------
figure(3);
localPlotRes(aTime,aSignal(:,2),aPar.aSample,'Geophone: 2, Cannel: Z');

%---------------------------------------------------------------------------
figure(4);
localPlotRes(aTime,aSignal(:,3),aPar.aSample,'Geophone: 3, Cannel: Z');

%---------------------------------------------------------------------------
figure(5);
Ad=80;
Ah=8;
AF1 = 0.5;
AF2 = 5.5;
for k=1:length(aTime)
    aDelta1(k)=asin((aSignal(k,1)-aSignal(k,2))/Ad);
    aDelta2(k)=asin((aSignal(k,2)-aSignal(k,3))/Ah);
end;
[aFreq3,aFourD1] = gwlFourTrans('mat',aDelta1,2,aSFreq);
[aFreq3,aFourD2] = gwlFourTrans('mat',aDelta2,2,aSFreq);
gwlPlotFunction(aFreq3,abs(aFourD1),0.07,0.2,0.4,0.4,AF1,AF2,0,max(abs(aFourD1)),gwlGetNotation('FREQ'),'F_\alpha','(a)');
gwlPlotFunction(aFreq3,abs(aFourD2),0.55,0.2,0.4,0.4,AF1,AF2,0,max(abs(aFourD2)),gwlGetNotation('FREQ'),'F_\beta','(b)');

%---------------------------------------------------------------------------
pause(0.00001);
print -f2 -r600 -depsc VibroBridgeFig2;
print -f3 -r600 -depsc VibroBridgeFig3;
print -f4 -r600 -depsc VibroBridgeFig4;
print -f5 -r600 -depsc VibroBridgeFig5;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function [aTime,aFreq,aWav] = localWavTrans(aSource,aFCount,aF1,aF2,aT1,aSampleFreq,aPar)
aTCount=length(aSource);
for n=1:aTCount
    aTime(n)=aT1 + (n-1)/aSampleFreq;
end;
aFour = fft(aSource);
for i=1:aFCount
    i
    aFreq(i)=aF1 + (i-1)*(aF2-aF1)/(aFCount-1);
    for n=1:aTCount
        afn = aSampleFreq*(n-1)/(aTCount-1);
        as0 = localWMorletF(afn/aFreq(i),aPar);
        aH(n) = aFour(n)*conj(as0);
    end;
    aSn = ifft(aH);
    for n=1:aTCount
        aWav(i,n) = sqrt(2*pi)*aSn(n);
    end;
end;


%---------------------------------------------------------------------------
function [aTime1,aFreq1,aWav1] = localWavCat(aTime,aFreq,aWav,aCols,aRows)
[aSizey,aSizex]=size(aWav);
aStepy = aSizey/aCols;
aStepx = aSizex/aRows;
m1=1;
k1=1;
for k=1:aRows
    aTime1(k) = aTime(k1);
    for m=1:aCols
        aWav1(m,k) = aWav(m1,k1);
        aFreq1(m) = aFreq(m1);
        m1 = m1+aStepy;
    end;
    k1 = k1+aStepx;
    m1 = 1;
end;


%---------------------------------------------------------------------------
function aWal = localWMorletF(aR, aPar)
aWal = aPar*exp(-(aPar*(2*pi*aR-2*pi)^2)/2);


%---------------------------------------------------------------------------
function localPlotRes(aTime1,aSignal,aSFreq,aNotation2)
aFCount = 256;
aYmax = 0.0002;
AT1 = min(aTime1);
AT2 = max(aTime1);
AF1 = 0.5;
AF2 = 5.5;
[aTime3,aFreq3,aWav3] = localWavTrans(aSignal,aFCount,AF1,AF2,AT1,aSFreq,1);
[aTime4,aFreq4,aWav4] = localWavCat(aTime3,aFreq3,aWav3,aFCount,256);
[aFreq,aFour] = gwlFourTrans('mat',aSignal,2,aSFreq);
gwlPlotFunction(aTime1,aSignal,0.07,0.75,0.9,0.20,AT1,AT2,-aYmax,aYmax,'',gwlGetNotation('SIG','T','z'),['(a) ',aNotation2]);
gwlPlotImage(aTime4,aFreq4,abs(aWav4),0.07,0.33,0.9,0.40,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('SIG','WABS','z')]);
gwlPlotFunction(aFreq,abs(aFour),0.32,0.07,0.39,0.20,AF1,AF2,0,max(abs(aFour)),gwlGetNotation('FREQ'),gwlGetNotation('SIG','FM','z'),'(c)');
