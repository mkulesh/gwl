function twoStations()
% twoStations(): Two component polarization analysis of a geomagnetic record.
% 
% In this example, we examine the phase differences of Pi2 pulsations in the 
% northward direction between HER (Hermanus, -33.9) and KAK stations for three 
% events occurred on 7 December 1996, 4 January 1997, and 14 January 1997. Figure 
% 1 shows geomagnetic field variations recorded at HER and KAK on the left panel, 
% and results of polarization analysis on the right panel. 
% 
% On the left panel for the 14 January 1997 event, we almost do not observe any 
% phase differences between two stations during the event. Pi2 pulsations appeared 
% in the time interval of 11.95-12.1h, during which oscillation peaks were 
% coincident between HER and KAK. This characteristic is reflected in the 
% polarization analysis shown on the right panel. There are not significant 
% changes in color, indicating no clear change in phase differences. The same 
% behavior we can observe also on the bottom panel of Figure 2, which corresponds 
% to this event. 
% 
% However, in the first two events (the 7 December 1996 event and 4 January 1997 
% event) we can see apparent phase change during events. At the beginning of 
% oscillations for the 1996-12-07 event (about 17.53h, region "A" on the top-right 
% panel in Figure 1), phase difference is about -50, that we can read from the 
% first panel of Figure 2. After 2-3 cycles of oscillations the phase difference 
% becomes positive value as indicated by the darker color in region "B". The 1997-
% 01-04 event demonstrates an opposed behavior: at the beginning of oscillations 
% the phase difference is slightly negative (region "C" on the middle-right panel 
% in Figure 1), but after 1-2 cycles it becomes negative value about -50 (region 
% "D") that can be identified by the change of color to more light. 
% 
% [1] M. Kulesh, M. Nose, M. Holschneider, K. Yumoto. Polarization analysis of a 
%     Pi2 pulsation using continuous wavelet transform // Preprint Series DFG SPP 
%     1114, University of Bremen. Preprint 153 (2007).
% [2] M. Kulesh, M. Nose and M. Holschneider. Polarization Analysis of Pi2 
%     Pulsations Using Continuous Wavelet Transform // Eos Trans. AGU, 87(52), Fall 
%     Meet. Suppl., Abstract SM43D-02 (2006).
% 
% FIGURE 1. North components of source signals and phase difference between HER 
% and KAK stations for three events: 1996-12-07, 1997-01-04 and 1997-01-14.
% 
% FIGURE 2. Time-dependent phase differencies averaged in the frequency band 0.01 
% and 0.015Hz.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aFreqName = 'freq.dat';
gwlCreateAxis(128,0.005,0.03,'lin --sign=full',aFreqName,'Frequency');
evalOneStation('1996-12-07', 'HER-KAK', '0,2', [62500 64548], aFreqName);
evalOneStation('1997-01-04', 'HER-KAK', '0,2', [56280 58328], aFreqName);
evalOneStation('1997-01-14', 'HER-KAK', '0,2', [42300 44348], aFreqName);
aGeom.aCat = 60;

%---------------------------------------------------------------------------
figure(1); 
aGeom.aText1 = 'A';
aGeom.aText1X = 17.53;
aGeom.aText1Y = 0.015;
aGeom.aText2 = 'B';
aGeom.aText2X = 17.65;
aGeom.aText2Y = 0.015;
plotsPolarizPhaseDiff('1996-12-07', 'HER-KAK', aGeom, 0,0.91,' ','HER','KAK',4);
aGeom.aText1 = 'C';
aGeom.aText1X = 15.76;
aGeom.aText1Y = 0.015;
aGeom.aText2 = 'D';
aGeom.aText2X = 15.86;
aGeom.aText2Y = 0.015;
plotsPolarizPhaseDiff('1997-01-04', 'HER-KAK', aGeom, 1,0.72,' ','HER','KAK',3);
aGeom.aText1 = '';
aGeom.aText2 = '';
plotsPolarizPhaseDiff('1997-01-14', 'HER-KAK', aGeom, 1,0.53,gwlGetNotation('TIME','hours'),'HER','KAK',5);


%---------------------------------------------------------------------------
figure(2); 
aGeom.aAver1=27;
aGeom.aAver2=52;
aEvents = cellstr(['1996-12-07';'1997-01-04';'1997-01-14']);
[aTime,aAver,aInfoNew] = calcAverPhaseDiff(aEvents,'HER-KAK',aGeom);
aPhaseMax = 90;
aInd1 = 1+aGeom.aCat;
aInd2 = length(aTime(:,1))-aGeom.aCat;
gwlPlotFunction(aTime(:,1),aAver(:,1),0.07,0.91-0.185,0.9,0.165,aTime(aInd1,1),aTime(aInd2,1),-aPhaseMax,aPhaseMax,' ',' ',aEvents(1));
gwlPlotFunction(aTime(:,2),aAver(:,2),0.07,0.72-0.185,0.9,0.165,aTime(aInd1,2),aTime(aInd2,2),-aPhaseMax,aPhaseMax,' ',gwlGetNotation('EPAR','PDIFF','T'),aEvents(2));
gwlPlotFunction(aTime(:,3),aAver(:,3),0.07,0.53-0.185,0.9,0.165,aTime(aInd1,3),aTime(aInd2,3),-aPhaseMax,aPhaseMax,gwlGetNotation('TIME','hours'),' ',aEvents(3));
aInfoNew

%---------------------------------------------------------------------------
pause(0.00001);
delete('*.dat');
clear all;
print -f1 -r600 -depsc twoStationsFig1;
print -f2 -r600 -depsc twoStationsFig2;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function plotsPolarizPhaseDiff(aFileName,aChanNot,aGeom,aColBar,aY,aXlab,aY1lab,aY2lab,aDeltaY)
aSignalName = strcat(aFileName,aChanNot,'sig.dat');
aElliparName = strcat(aFileName,aChanNot,'par.dat');
fid = fopen(aSignalName,'r');    [aTime,aSig] = gwlReadSignal(fid);    fclose(fid);  
fid = fopen(aElliparName,'r');   [aTime,aFreq,aCWT] = gwlReadSpectrum(fid);   fclose(fid);  
aTime = aTime./3600;
aInd1 = 1+aGeom.aCat;
aInd2 = length(aTime)-aGeom.aCat;
aPhDiff = aCWT(:,:,1);
aPhDiff(1,aInd1+1) = 180;
aPhDiff(1,aInd1+2) = -180;
aX = 0.05;
aTitle = '';
if(aColBar == 0)
    aTitle = gwlGetNotation('CSIG','T');
end;
aSig1 = real(aSig)-min(real(aSig));
aSig2 = imag(aSig)-min(imag(aSig));

gwlPlotFunction(aTime,aSig1,aX,aY-0.185,0.43,0.165,aTime(aInd1),aTime(aInd2),0,aDeltaY,aXlab,' ',aFileName,aTitle);
    hold on; plot(aTime,aSig2,'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend(aY1lab,aY2lab);
if(aColBar == 0)
    aLength = 0.493;
else  
    aLength = 0.42;
end;  
if(aColBar == 0)
    aTitle = gwlGetNotation('EPAR','PDIFF','WG');
end;
gwlPlotImage(aTime(aInd1:aInd2),aFreq,aPhDiff(:,aInd1:aInd2),aX+0.44,aY-0.185,aLength,0.165,aXlab,gwlGetNotation('FREQ'),'',aTitle);
    set(gca,'Box','Off');
    set(gca,'YAxisLocation','right');
    if(aColBar == 0)
        set(gca,'YTickLabel',{});
        ylabel('');
        clb = colorbar;
        set(clb,'FontSize',7);
    end;  
    gwlText(aGeom.aText1X,aGeom.aText1Y,aGeom.aText1);
    gwlText(aGeom.aText2X,aGeom.aText2Y,aGeom.aText2);

%---------------------------------------------------------------------------
function [aTime,aAver,aInfoNew] = calcAverPhaseDiff(aEvents,aStation,aGeom)
for n=1:length(aEvents)
    aElliparName = strcat(char(aEvents(n)),aStation,'par.dat');
    fid = fopen(aElliparName,'r');   [aTimeTmp,aFreq,aCWT] = gwlReadSpectrum(fid);   fclose(fid);  
    aTime(:,n) = aTimeTmp(:,1)./3600;
    aPhDiff = aCWT(:,:,1);
    for m=1:length(aTime)
        aSumPoint=1;
        aVal = 0;
        for k=aGeom.aAver1:aGeom.aAver2
            if aPhDiff(k,m)~=0
                aVal = aVal+aPhDiff(k,m);
                aSumPoint = aSumPoint+1;
            end
        end
        aAver(m,n) = aVal/aSumPoint;
    end
end
aInfoNew.aAver1 = aFreq(aGeom.aAver1);
aInfoNew.aAver2 = aFreq(aGeom.aAver2);

