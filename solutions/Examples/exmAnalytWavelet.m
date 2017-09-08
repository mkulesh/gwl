function exmAnalytWavelet()


%---------------------------------------------------------------------------
path(path, '../../mshell');
aPar.aSurfCount = 50;
[aTime,aPart] = gwlCreateAxis(600,0,1);
aPar.aSample = aPart.aSample;
aFreq = gwlCreateAxis(300,0.001,50);
aFreqDouble = gwlCreateAxis(300,0.001,100);


%---------------------------------------------------------------------------
figure(1)
aT0 = aTime(length(aTime)/2);
for n=1:length(aTime)
    aSig(n) = localSamp1(0,aTime(n),0,0.2,aT0);
    for m=1:length(aFreq)
        aWav(m,n) = localSamp1(1,aTime(n),aFreq(m),0.2,aT0);
    end;
end;
localPlotRes(aTime, aFreq, aSig, aWav, aPar, gwlGetNotation('SIG','T',1), gwlGetNotation('SIG','FM',1), 0.00005, 'NoFour');
 
%---------------------------------------------------------------------------
figure(2)
for n=1:length(aTime)
    aSig(n) = localSamp2(0,aTime(n),0,2,25);
    for m=1:length(aFreq)
        aWav(m,n) = localSamp2(1,aTime(n),aFreq(m),2,25);
    end;
end;
localPlotRes(aTime, aFreq, aSig, aWav, aPar, gwlGetNotation('SIG','T',2), gwlGetNotation('SIG','FM',2));

%---------------------------------------------------------------------------
figure(3)
for n=1:length(aTime)
    aSig(n) = localSamp3(0,aTime(n),0,3,15,10);
    for m=1:length(aFreq)
        aWav(m,n) = localSamp3(1,aTime(n),aFreq(m),3,15,10);
    end;
end;
localPlotRes(aTime, aFreq, aSig, aWav, aPar, gwlGetNotation('SIG','T',3), gwlGetNotation('SIG','FM',3));

%---------------------------------------------------------------------------
figure(4)
for n=1:length(aTime)
    aSig(n) = localSamp3(0,aTime(n),0,3,15,10) + localSamp3(0,aTime(n),0,3,35,-10);
    for m=1:length(aFreq)
        aWav(m,n) = localSamp3(1,aTime(n),aFreq(m),3,15,10) + localSamp3(1,aTime(n),aFreq(m),3,35,-10);
    end;
end;
localPlotRes(aTime, aFreq, aSig, aWav, aPar, gwlGetNotation('SIG','T',4), gwlGetNotation('SIG','FM',4));

%---------------------------------------------------------------------------
figure(5)
for n=1:length(aTime)
    aSig(n) = localSamp4(0,aTime(n),0,1,0.5,25,1);
    for m=1:length(aFreq)
        aWav(m,n) = localSamp4(1,aTime(n),aFreq(m),1,0.5,25,1);
    end;
end;
localPlotRes(aTime, aFreq, aSig, aWav, aPar, gwlGetNotation('SIG','T',5), gwlGetNotation('SIG','FM',5));

%---------------------------------------------------------------------------
figure(6)
for n=1:length(aTime)
    aSig(n) = localSamp4(0,aTime(n),0,1,0.25,10,1) + localSamp4(0,aTime(n),0,1,0.6,30,1) + localSamp4(0,aTime(n),0,1,0.8,70,1);
    for m=1:length(aFreqDouble)
        aWav(m,n) = localSamp4(1,aTime(n),aFreqDouble(m),1,0.25,10,1) + localSamp4(1,aTime(n),aFreqDouble(m),1,0.6,30,1) + localSamp4(1,aTime(n),aFreqDouble(m),1,0.8,70,1);
    end;
end;
localPlotRes(aTime, aFreqDouble, aSig, aWav, aPar, gwlGetNotation('SIG','T',6), gwlGetNotation('SIG','FM',6));
 
%---------------------------------------------------------------------------
figure(7)
for n=1:length(aTime)
    aSig(n) = localSamp4(0,aTime(n),0,1,0.25,10,1)/10 + localSamp4(0,aTime(n),0,1,0.6,30,1)/30 + localSamp4(0,aTime(n),0,1,0.8,70,1)/70;
    for m=1:length(aFreqDouble)
        aWav(m,n) = localSamp4(1,aTime(n),aFreqDouble(m),1,0.25,10,1)/10 + localSamp4(1,aTime(n),aFreqDouble(m),1,0.6,30,1)/30 + localSamp4(1,aTime(n),aFreqDouble(m),1,0.8,70,1)/70;
    end;
end;
localPlotRes(aTime, aFreqDouble, aSig, aWav, aPar, gwlGetNotation('SIG','T',7), gwlGetNotation('SIG','FM',7));

%---------------------------------------------------------------------------
figure(8)
for n=1:length(aTime)
    aSig(n) = localSamp5(0,aTime(n),0,1,0.5,50);
    for m=1:length(aFreq)
      aWav(m,n) = localSamp5(1,aTime(n),aFreq(m),1,0.5,50);
  end;
end;
localPlotRes(aTime, aFreq, aSig, aWav, aPar, gwlGetNotation('SIG','T',8), gwlGetNotation('SIG','FM',8));

%---------------------------------------------------------------------------
pause(0.00001);
clear all;
print -f1 -r600 -depsc exmAnalytWaveletFig1;
print -f2 -r600 -depsc exmAnalytWaveletFig2;
print -f3 -r600 -depsc exmAnalytWaveletFig3;
print -f4 -r600 -depsc exmAnalytWaveletFig4;
print -f5 -r600 -depsc exmAnalytWaveletFig5;
print -f6 -r600 -depsc exmAnalytWaveletFig6;
print -f7 -r600 -depsc exmAnalytWaveletFig7;
print -f8 -r600 -depsc exmAnalytWaveletFig8;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function  localPlotRes(aR, aF, aS, aW, aPar, aNotation, aNotationFreq, aEps, aNoFour)
if (nargin < 8) 
    aEps = 0.1;
end;
[aFreq1,aFour1] = gwlFourTrans('mat',aS,2,aPar.aSample);
aYmax = max(real(aS))*1.2;
gwlPlotFunction(aR,real(aS),0.07,0.75,0.9,0.20,min(aR),max(aR),-aYmax,aYmax,'',aNotation,'(a)');
    hold on; plot(aR,imag(aS),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
gwlPlotSurface(aR,aF,aW,0,0,aPar.aSurfCount,0.07,0.54,0.9,0.20,'',gwlGetNotation('FREQ'),'(b)');
    set(gca,'Box','On');
gwlPlotSurface(aR,aF,aW,1,aEps,aPar.aSurfCount,0.07,0.33,0.9,0.20,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(c)');
    set(gca,'Box','On');
if (nargin < 9) 
    gwlPlotFunction(aFreq1,abs(aFour1),0.35,0.07,0.33,0.20,0.0,max(aF),0,max(abs(aFour1)),gwlGetNotation('FREQ'),aNotationFreq,'(d)');
end;    

%---------------------------------------------------------------------------
function  aWal = localSamp1(aType, aB, aF, aPar, aTau)
if(aType == 0)
    if(aB == aTau)  
        aWal = 1.0;
    else
        aWal = 0.0;
    end;
else
    aMod = aF*exp(-(aF*(aTau-aB)/aPar)^2/2.0);
    aPhase = -aF*(aTau-aB);
    aWal = aMod*exp(2*pi*i*aPhase);
end;

%---------------------------------------------------------------------------
function  aWal = localSamp2(aType, aB, aF, aPar, aOm)
if(aType == 0)
    aWal = exp(2*pi*i*aOm*aB);
else
    aMod = aPar*sqrt(2*pi)*exp(-2*(pi*aPar*(aF-aOm)/aF)^2);
    aPhase = aB*aOm;
    aWal = aMod*exp(2*pi*i*aPhase);
end;

%---------------------------------------------------------------------------
function  aWal = localSamp3(aType, aB, aF, aPar, aOm, aDom)
if(aType == 0)
    aWal = exp(2*pi*i*(aOm+aDom*aB)*aB);
else
    aMod = 2*aF*aPar*sqrt(pi)/sqrt(2*aF^2 - 8*pi*i*aPar^2*aDom);
    aPhase = (-aDom*(aF*aB)^2 + 4*pi*i*aF*aB*aDom*aPar^2 - pi*i*(aOm*aPar)^2 - aOm*aB*aF^2 + 2*pi*i*aOm*aF*aPar^2 - pi*i*(aF*aPar)^2)/(-aF^2 + 4*pi*i*aDom*aPar^2);
    aWal = aMod*exp(2*pi*i*aPhase);
end;

%---------------------------------------------------------------------------
function  aWal = localSamp4(aType, aB, aF, aPar, aBm, aFm, aParm)
if(aType == 0)
    aWal = aFm*exp(-((aFm*(aB-aBm)/aParm)^2)/2.0)*exp(2*pi*i*aFm*(aB-aBm));
else
    aMod1 = exp(-((aF*aB*aFm)^2 + (aFm*aBm*aF)^2 + (2*pi*aF*aPar*aParm)^2 + (2*pi*aFm*aPar*aParm)^2 - 2*aB*aBm*(aF*aFm)^2 - 8*aF*aFm*(pi*aPar*aParm)^2)/(2*(aF*aParm)^2 + 2*(aFm*aPar)^2));
    aMod2 = aPar*aParm*aF*aFm*sqrt(2*pi)/sqrt((aF*aParm)^2 + (aFm*aPar)^2);
    aPhase = -aF*aFm*(aF*aBm*aParm^2 - aB*aFm*aPar^2 + aFm*aBm*aPar^2 - aF*aB*aParm^2)/((aF*aParm)^2 + (aFm*aPar)^2);
    aWal = aMod1*aMod2*exp(2*pi*i*aPhase);
end;

%---------------------------------------------------------------------------
function  aWal = localSamp5(aType, aB, aF, aPar, aBm, aFm)
if(aType == 0)
    aWal = aFm*(1-((aB-aBm)*aFm)^2)*exp(-(((aB-aBm)*aFm)^2)/2);
else
    aMod1 = exp(-((aF^2)/2.0)*((aB*aFm)^2 + (aFm*aBm)^2 + (2*pi*aPar)^2 - 2*aB*aBm*(aFm^2))/((aFm*aPar)^2 + aF^2));
    aMod2 = 8*i*sqrt(2*pi)*pi*(aPar^3)*(aF^4)*(aFm^3)*(aB-aBm)/sqrt(((aFm*aPar)^2 + aF^2)^5);
    aMod3 = 2*sqrt(2*pi)*aPar*(aF^3)*aFm*((aFm*aPar)^2 + aF^2 + (2*pi*aFm*aPar)^2 - (aFm*aF*aB)^2 + 2*aBm*aB*(aF*aFm)^2 - (aFm*aBm*aF)^2)/sqrt(((aFm*aPar)^2 + aF^2)^5);
    aPhase = aF*((aPar*aFm)^2)*(aB-aBm)/((aFm*aPar)^2 + aF^2);
    aWal = aMod1*(aMod2+aMod3)*exp(2*pi*i*aPhase);
end;

