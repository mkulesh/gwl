function SynthSig() 
% In this example, we test the adaptive nearest neighbours method (gwlNNpred) 
% on three deterministic time series with pronounced season dependence, trend 
% and amplitude change. Using adaptive  cityblock metrics, we calculate the 
% prediction and the MAPE value for each time series. For this calculation, 
% we select the embedding dimension similar to the seasonal period. 
%
% [1] K. Kurennaya, M. Kulesh, M. Holschneider Adaptive metrics in the nearest 
%     neighbours method // Preprint Series DFG SPP 1114, University of Bremen. 
%     Preprint 139 (2006).
% [2] M. Kulesh, et al., Adaptive metrics in the nearest neighbours 
%     method, Physica D (2007), doi:10.1016/j.physd.2007.08.019
%
% FIGURE 1. (a) Time series with season dependence, vertical line shows the 
% prediction start. (b) Zoom of predicted values and (c) Fourier spectrum of 
% source data
%
% FIGURE 2. (a) Time series with multiplicative seasonality, vertical line 
% shows the prediction start. (b) Zoom of predicted values and (c) Fourier 
% spectrum of source data
%
% FIGURE 3. (a) High-frequency time series with multiplicative seasonality, 
% vertical line shows the prediction start. (b) Zoom of predicted values and 
% (c) Fourier spectrum of source data

%----------------------------------------------------------------------------
path(path, '../../mshell');

%----------------------------------------------------------------------------
figure(1); 
localPlotOneSerie('SynthSig1.asc', 'SynthSig1.dat',0.02,0); 

%----------------------------------------------------------------------------
figure(2); 
localPlotOneSerie('SynthSig2.asc', 'SynthSig2.dat',0.5,359); 

%----------------------------------------------------------------------------
figure(3); 
localPlotOneSerie('SynthSig3.asc', 'SynthSig3.dat',0.5,0); 

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc SynthSigFig1;
print -f2 -r600 -depsc SynthSigFig2;
print -f3 -r600 -depsc SynthSigFig3;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function localPlotOneSerie(aSourceFile, aResFile, aFreqMax, aT0) 
aSourceFile
aSource = load(aSourceFile,'-ascii');
aRes    = load(aResFile,'-ascii'); 
aStart  = length(aSource)-length(aRes)+1;
aEnd    = length(aSource);
calcMAPE(aSource, aRes, aStart, aEnd)

aMin = 0;
aMax = max(aSource);
aAmpl=aMin:(aMax-aMin)/100:aMax;
aTime(1,:) = (1:aEnd)+aT0;
aTime1(1,:) = (aStart:aEnd)+aT0;
gwlPlotFunction(aTime,aSource,0.07,0.57,0.9,0.4,min(aTime),max(aTime),aMin,aMax,'Time','Data series','(a)');
    hold on; 
    plot(aTime1,aRes,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);
    hold off;
    line([aStart+aT0,aStart+aT0],[aMin,aMax],'Color',gwlGetColor(1));
    grid off;

aMin = min(aSource(aStart:aEnd));
aMax = max(aSource(aStart:aEnd));
gwlPlotFunction(aTime1,aSource(aStart:aEnd),0.07,0.07,0.4,0.4,min(aTime1),max(aTime1),aMin,aMax*1.1,'Time','Data series','(b)');
    hold on; 
    plot(aTime1,aRes,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);
    hold off;
    legend('Source data','Prediction using L_A');
    grid off;
[aFreq,aFour] = gwlFourTrans('mat',aSource,2,1);
aFourAbs = log(abs(aFour));
gwlPlotFunction(aFreq,aFourAbs,0.57,0.07,0.4,0.4,0,aFreqMax,0,max(aFourAbs),'Frequency','log(Power)','(c)');
    grid off;

