function eventSeptember20()
% eventSeptember20(): Two component polarization analysis of a geomagnetic record.
%
% At 0538 UT on 20 September 1995, Pi2 pulsations were observed clearly at three 
% different MLT sectors. At postmidnight sector (0200 MLT), the low-latitude 
% station SMA (Santa Maria, -19.5 geomagnetic latitude (GMLAT)) and the 
% equatorial station BLM (Belem, 8.6 GMLAT) observed clear Pi2 pulsations. On the 
% dawnside (0700 MLT), Pi2 pulsation was detected at the mid-latitude station LAQ 
% (Laquila, 42.5 GMLAT). Even in the afternoon sector at 1500 MLT, Pi2 
% pulsations were found at low-latitude stations MSR (Moshiri, 35.4 GMLAT) and 
% KAK (Kakioka, 26.9 GMLAT), and at the equatorial station GUA (Guam, 5.1 
% GMLAT). This event was reported by Nose et. al, Multipoint observations of a Pi2 
% pulsation on morningside: The 20 September 1995 event. Journal of Geophysical
% Research 108 (A5), 1219 (2003). They suggested that the substorm was initiated 
% at pre-midnight (2100-2300 MLT) by using high-latitude ground magnetic field 
% data. Geomagnetic field data in the northward and eastward components recorded 
% at MSR, GUA, SMA, BLM, and LAQ is demonstrated on the Figure 1.
% 
% [1] M. Kulesh, M. Nose, M. Holschneider, K. Yumoto. Polarization analysis of a 
%     Pi2 pulsation using continuous wavelet transform // Preprint Series DFG SPP 
%     1114, University of Bremen. Preprint 153 (2007).
% [2] M. Kulesh, M. Nose and M. Holschneider. Polarization Analysis of Pi2 
%     Pulsations Using Continuous Wavelet Transform // Eos Trans. AGU, 87(52), Fall 
%     Meet. Suppl., Abstract SM43D-02 (2006).
%     
% FIGURE 1. Source signals for different stations.
% 
% FIGURE 2. Source signals and time-frequency representations of phase 
% differences. Vertical lines bound the interval of Pi2 pulsation.    
% 
% FIGURE 3. Ellipticity ratio and tilt angle in time-frequency domain for signals 
% shown in Figure 2.
% 
% FIGURE 4 and FIGURE 5. Northward components of source signals registered at 
% different stations and time-frequency phase difference between the pair of 
% stations.
% 
% FIGURE 6. Frequency-dependent phase differencies averaged in the time interval 
% between 5:64h and 5:66h.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aFreqName = 'freq.dat';
gwlCreateAxis(128,0.005,0.06,'lin --sign=full',aFreqName,'Frequency');
evalOneStation('1995-09-20','MSR',     '2,3',  [19800 21336], aFreqName);
evalOneStation('1995-09-20','GUA',     '6,7',  [19800 21336], aFreqName);
evalOneStation('1995-09-20','SMA',     '8,9',  [19800 21336], aFreqName);
evalOneStation('1995-09-20','BLM',     '10,11',[19800 21336], aFreqName);
evalOneStation('1995-09-20','LAQ',     '12,13',[19800 21336], aFreqName);

aGeom.aInd1 = 80;
aGeom.aInd2 = 300;
aGeom.aPiLeft = 150;
aGeom.aPiRight = 250;
aGeom.aAver1 = 175;
aGeom.aAver2 = 195;
aGeom.aValTime = 180;
aGeom.aValFreq = 32;

%---------------------------------------------------------------------------
figure(1); 
SourceSig = load('1995-09-20.asc','-ascii');
aTime = SourceSig(:,1)/3600;
aDeltaY = 10;
aYmin = [-5 -54 -5 5 -2 32];
aInd = [2 4 6 8 10 12];
aText = cellstr(['HER, Morning  ';'MSR, Afternoon';'KAK, Afternoon';'GUA, Afternoon';'SMA, Midnight ';'BLM, Midnight ']);
aButtom = [0.90 0.81 0.72 0.63 0.54 0.45];
for k=1:6
    aXTitle = '';
    if k==6
        aXTitle = gwlGetNotation('TIME','hours');
    end;
    gwlPlotFunction(aTime,SourceSig(:,aInd(k)),0.2,aButtom(k),0.6,0.08,min(aTime),max(aTime),aYmin(k),aYmin(k)+aDeltaY,aXTitle,' ',aText(k));
    set(gca,'YTick',aYmin(k):(aDeltaY/5):(aYmin(k)+aDeltaY));
    set(gca,'YTickLabel',{' ','2nT','4nT','6nT','8nT',' '});
end;

%---------------------------------------------------------------------------
figure(2); 
aInfo = plotsPolarizPhaseDiff('1995-09-20','MSR', aGeom, 0,0.980,'','MSR, North','MSR, East',-54,33);
plotsPolarizPhaseDiff('1995-09-20','GUA', aGeom, 1,0.795,'','GUA, North','GUA, East',6,-34);
plotsPolarizPhaseDiff('1995-09-20','SMA', aGeom, 1,0.610,'','SMA, North','SMA, East',-1,-10);
plotsPolarizPhaseDiff('1995-09-20','BLM', aGeom, 1,0.425,'','BLM, North','BLM, East',32,6);
plotsPolarizPhaseDiff('1995-09-20','LAQ', aGeom, 1,0.240,gwlGetNotation('TIME','hours'),'LAQ, North','LAQ, East',0,2);

%---------------------------------------------------------------------------
figure(3); 
plotsPolarizRatioTilt('1995-09-20','MSR', aGeom, 0,0.980,'','MSR, North','MSR, East');
plotsPolarizRatioTilt('1995-09-20','GUA', aGeom, 1,0.795,'','GUA, North','GUA, East');
plotsPolarizRatioTilt('1995-09-20','SMA', aGeom, 1,0.610,'','SMA, North','SMA, East');
plotsPolarizRatioTilt('1995-09-20','BLM', aGeom, 1,0.425,'','BLM, North','BLM, East');
plotsPolarizRatioTilt('1995-09-20','LAQ', aGeom, 1,0.240,gwlGetNotation('TIME','hours'),'LAQ, North','LAQ, East');

%---------------------------------------------------------------------------
figure(4); 
evalOneStation('1995-09-20','LAQ-MSR', '12,2', [19800 21336], aFreqName);
evalOneStation('1995-09-20','LAQ-GUA', '12,6', [19800 21336], aFreqName);
evalOneStation('1995-09-20','LAQ-SMA', '12,8', [19800 21336], aFreqName);
evalOneStation('1995-09-20','LAQ-BLM', '12,10',[19800 21336], aFreqName);
plotsPolarizPhaseDiff('1995-09-20','LAQ-MSR', aGeom, 0,0.980,'','LAQ','MSR',0,-54);
plotsPolarizPhaseDiff('1995-09-20','LAQ-GUA', aGeom, 1,0.795,'','LAQ','GUA',0,6);
plotsPolarizPhaseDiff('1995-09-20','LAQ-SMA', aGeom, 1,0.610,'','LAQ','SMA',0,-1);
plotsPolarizPhaseDiff('1995-09-20','LAQ-BLM', aGeom, 1,0.425,'','LAQ','BLM',0,32);

%---------------------------------------------------------------------------
figure(5); 
evalOneStation('1995-09-20','MSR-GUA', '2,6', [19800 21336], aFreqName);
evalOneStation('1995-09-20','MSR-SMA', '2,8', [19800 21336], aFreqName);
evalOneStation('1995-09-20','MSR-BLM', '2,10',[19800 21336], aFreqName);
plotsPolarizPhaseDiff('1995-09-20','MSR-GUA', aGeom, 1,0.795,'','MSR','GUA',-54,4);
plotsPolarizPhaseDiff('1995-09-20','MSR-SMA', aGeom, 1,0.610,'','MSR','SMA',-54,-4);
plotsPolarizPhaseDiff('1995-09-20','MSR-BLM', aGeom, 1,0.425,gwlGetNotation('TIME','hours'),'MSR','BLM',-54,29);

%---------------------------------------------------------------------------
figure(6); 
aStation = cellstr(['MSR';'GUA';'SMA';'BLM';'LAQ']);
[aFreq,aAver,aInfoNew] = calcAverPhaseDiff('1995-09-20',aStation,aGeom,aInfo);
gwlPlotFunction(aFreq,aAver(:,1),0.07,0.07,0.9,0.9,min(aFreq),max(aFreq),-180,180,gwlGetNotation('FREQ'),gwlGetNotation('EPAR','PDIFF','F'));
    hold on;  
    aTxtIndex = 30;
    plot(aFreq,aAver(:,2),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);
    plot(aFreq,aAver(:,3),'Color',gwlGetColor(0),'LineStyle',':','LineWidth',2);
    plot(aFreq,aAver(:,4),'Color',gwlGetColor(0),'LineStyle','-.','LineWidth',1);
    plot(aFreq,aAver(:,5),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    hold off;
    for n=1:length(aStation)
        gwlText(aFreq(aTxtIndex),aAver(aTxtIndex,n),char(aStation(n)));
    end;

%---------------------------------------------------------------------------
pause(0.00001);
aInfoNew
delete('*.dat');
clear all;

print -f1 -r600 -depsc eventSeptember20Fig1;
print -f2 -r600 -depsc eventSeptember20Fig2;
print -f3 -r600 -depsc eventSeptember20Fig3;
print -f4 -r600 -depsc eventSeptember20Fig4;
print -f5 -r600 -depsc eventSeptember20Fig5;
print -f6 -r600 -depsc eventSeptember20Fig6;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function aInfo = plotsPolarizPhaseDiff(aFileName,aChanNot,aGeom,aColBar,aY,aXlab,aY1lab,aY2lab,aYmin1,aYmin2)
aSignalName = strcat(aFileName,aChanNot,'sig.dat');
aElliparName = strcat(aFileName,aChanNot,'par.dat');
fid = fopen(aSignalName,'r');    [aTime,aSig] = gwlReadSignal(fid);    fclose(fid);  
fid = fopen(aElliparName,'r');   [aTime,aFreq,aCWT] = gwlReadSpectrum(fid);   fclose(fid);  
aTime = aTime./3600;
aPhDiff = aCWT(:,:,1);
aPhDiff(1,aGeom.aInd1+1) = 180;
aPhDiff(1,aGeom.aInd1+2) = -180;
aX = 0.05;
aDeltaY = 10;
aTitle = '';
if(aColBar == 0)
    aTitle = gwlGetNotation('CSIG','T');
end;
gwlPlotFunction(aTime,real(aSig)-aYmin1,aX,aY-0.185,0.43,0.165,aTime(aGeom.aInd1),aTime(aGeom.aInd2),0,aDeltaY,aXlab,' ','',aTitle);
   hold on; plot(aTime,imag(aSig)-aYmin2,'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
   line([aTime(aGeom.aPiLeft),aTime(aGeom.aPiLeft)],[0,(aDeltaY)],'Color',gwlGetColor(2));
   line([aTime(aGeom.aPiRight),aTime(aGeom.aPiRight)],[0,(aDeltaY)],'Color',gwlGetColor(2));
   set(gca,'YTick',0:(aDeltaY/5):(aDeltaY));
   set(gca,'YTickLabel',{'0','2nT','4nT','6nT','8nT',' '});
   legend(aY1lab,aY2lab);
   
   if(aColBar == 0)
       aLength = 0.493;
   else  
       aLength = 0.42;
   end;  
aTitle = '';
if(aColBar == 0)
    aTitle = gwlGetNotation('EPAR','PDIFF','WG');
end;
gwlPlotImage(aTime(aGeom.aInd1:aGeom.aInd2),aFreq,aPhDiff(:,aGeom.aInd1:aGeom.aInd2),aX+0.44,aY-0.185,aLength,0.165,aXlab,gwlGetNotation('FREQ'),aChanNot,aTitle);
    line([aTime(aGeom.aPiLeft),aTime(aGeom.aPiLeft)],[min(aFreq),max(aFreq)],'Color','white');
    line([aTime(aGeom.aPiRight),aTime(aGeom.aPiRight)],[min(aFreq),max(aFreq)],'Color','white');
    set(gca,'Box','Off');
    set(gca,'YAxisLocation','right');
    if(aColBar == 0)
        set(gca,'YTickLabel',{});
        ylabel('');
        clb = colorbar;
        set(clb,'FontSize',7);
    end;  
    line([aTime(aGeom.aValTime),aTime(aGeom.aValTime)*1.013],[aFreq(aGeom.aValFreq),aFreq(aGeom.aValFreq)*1.7],'Color','white');
    gwlText(aTime(aGeom.aValTime)*1.012,aFreq(aGeom.aValFreq)*2.1,num2str(aPhDiff(aGeom.aValFreq,aGeom.aValTime)));
aInfo.aValTime = aTime(aGeom.aValTime);
aInfo.aValFreq = aFreq(aGeom.aValFreq);

%---------------------------------------------------------------------------
function plotsPolarizRatioTilt(aFileName,aChanNot,aGeom,aColBar,aY,aXlab,aY1lab,aY2lab)
aSignalName = strcat(aFileName,aChanNot,'sig.dat');
aElliparName = strcat(aFileName,aChanNot,'par.dat');
fid = fopen(aSignalName,'r');    [aTime,aSig] = gwlReadSignal(fid);    fclose(fid);  
fid = fopen(aElliparName,'r');   [aTime,aFreq,aCWT] = gwlReadSpectrum(fid);   fclose(fid);  
aTime = aTime./3600;
aRatio = aCWT(:,:,2);
aRatio(1,aGeom.aInd1+1) = 0;
aRatio(1,aGeom.aInd1+2) = 1;
aTilt = aCWT(:,:,3);
aTilt(1,aGeom.aInd1+1) = 90;
aTilt(1,aGeom.aInd1+2) = -90;
aX = 0.01;
if(aColBar == 0)
    aLength = 0.463;
else  
    aLength = 0.395;
end;  
aTitle = '';
if(aColBar == 0)
    aTitle = gwlGetNotation('EPAR','RATIO','WG');
end;
gwlPlotImage(aTime(aGeom.aInd1:aGeom.aInd2),aFreq,aRatio(:,aGeom.aInd1:aGeom.aInd2),aX,aY-0.185,aLength,0.165,aXlab,gwlGetNotation('FREQ'),aChanNot,aTitle);
    line([aTime(aGeom.aPiLeft),aTime(aGeom.aPiLeft)],[min(aFreq),max(aFreq)],'Color',gwlGetColor(2));
    line([aTime(aGeom.aPiRight),aTime(aGeom.aPiRight)],[min(aFreq),max(aFreq)],'Color',gwlGetColor(2));
    set(gca,'Box','Off');
    set(gca,'YAxisLocation','right');
    if(aColBar == 0)
        set(gca,'YTickLabel',{});
        ylabel('');
        clb = colorbar;
        set(clb,'FontSize',7);
    end;  
aTitle = '';
if(aColBar == 0)
    aTitle = gwlGetNotation('EPAR','TILT','WG');
end;
gwlPlotImage(aTime(aGeom.aInd1:aGeom.aInd2),aFreq,aTilt(:,aGeom.aInd1:aGeom.aInd2),aX+0.48,aY-0.185,aLength,0.165,aXlab,gwlGetNotation('FREQ'),aChanNot,aTitle);
    line([aTime(aGeom.aPiLeft),aTime(aGeom.aPiLeft)],[min(aFreq),max(aFreq)],'Color','white');
    line([aTime(aGeom.aPiRight),aTime(aGeom.aPiRight)],[min(aFreq),max(aFreq)],'Color','white');
    set(gca,'Box','Off');
    set(gca,'YAxisLocation','right');
    if(aColBar == 0)
        set(gca,'YTickLabel',{});
        ylabel('');
        clb = colorbar;
        set(clb,'FontSize',7);
    end;  
   

%---------------------------------------------------------------------------
function [aFreq,aAver,aInfoNew] = calcAverPhaseDiff(aFileName,aStation,aGeom,aInfo)
for n=1:length(aStation)
    aElliparName = strcat('1995-09-20',char(aStation(n)),'par.dat');
    fid = fopen(aElliparName,'r');   [aTime,aFreq,aCWT] = gwlReadSpectrum(fid);   fclose(fid);  
    aPhDiff = aCWT(:,:,1);
    for m=1:length(aFreq)
        aSumPoint=1;
        aVal = 0;
        for k=aGeom.aAver1:aGeom.aAver2
            if aPhDiff(m,k)~=0
                aVal = aVal+aPhDiff(m,k);
                aSumPoint = aSumPoint+1;
            end
        end
        aAver(m,n) = aVal/aSumPoint;
    end
end
aInfoNew = aInfo;
aInfoNew.aAver1 = aTime(aGeom.aAver1)./3600;
aInfoNew.aAver2 = aTime(aGeom.aAver2)./3600;