function filterHER()
% filterHER(): Two component polarization filter of a geomagnetic record.
%
% In this example, we show an example of polarization filter calculated using the 
% module gwlET2DFilter. We consider here three different events showed above in 
% the Fig. 1, but now we analyze the whole signal within the day containing a Pi2 
% pulsation. In Fig. 2, with the aim to extract the Pi2 pulsation from all these 
% signals, we use Fourier band-pass filter and the so-called "total horizontal 
% filter", which is defined in the same frequency band as Fourier filter, but with 
% additional restriction on the tilt angle |\tilt|< 0.6 rad. We can see that the 
% polarization filter picks out the Pi2 pulsation clearer than Fourier band-pass 
% filter. The next example in Fig. 3 is a division of previous horizontal 
% polarized signals into the linear and elliptical parts, where \ratio_f=0.3. As 
% expected, all of Pi2 pulsations are placed only on the left panels demonstrated 
% linearly polarized signals, because Pi2 pulsations on the nightside are 
% predominantly linearly polarized.
% 
% [1] M. Kulesh, M. Nose and M. Holschneider. Polarization Analysis of Pi2 
%     Pulsations Using Continuous Wavelet Transform // Eos Trans. AGU, 87(52), Fall 
%     Meet. Suppl., Abstract SM43D-02 (2006).
% 
% [2] M. Kulesh, M. Holschneider and M.S. Diallo. Geophysics Wavelet Library: 
%     Applications of the Continuous Wavelet Transform to the Polarization and 
%     Dispersion Analysis of Signals. Preprint series of the DFG priority program 1114 
%     "Mathematical methods for time series analysis and digital image processing", 
%     Preprint 156 (2007). 
% 
% FIGURE 1. Polarization properties of three geomagnetic records.
% 
% FIGURE 2. An example of total horizontal filter as applied to geomagnetic records.
% 
% FIGURE 3. An example of linear and elliptical filters as applied to geomagnetic records.

%---------------------------------------------------------------------------
path(path, '../../mshell');

aFreqName = 'freq.dat';
gwlCreateAxis(128,0.009,0.025,'lin --sign=full',aFreqName,'Frequency');
aGeom.aCat = 40;

evalOneStation('1996-12-07', 'HER', '0,1', [62500 64548], aFreqName);
evalOneStation('1997-01-04', 'HER', '0,1', [56280 58328], aFreqName);
evalOneStation('1997-01-14', 'HER', '0,1', [42300 44348], aFreqName);
[aTime,aSigBand1,aSigPolarLin1,aSigPolarElli1] = evalFilterRes(aFreqName,'1996-12-07', '0,1', 'HER', [10000 75536],'linhor,1.6,1,linhor,0.7,0.3,ellihor,0.7,0.3');
[aTime,aSigBand2,aSigPolarLin2,aSigPolarElli2] = evalFilterRes(aFreqName,'1997-01-04', '0,1', 'HER', [10000 75536],'linhor,1.6,1,linhor,0.6,0.3,ellihor,0.6,0.3');
[aTime,aSigBand3,aSigPolarLin3,aSigPolarElli3] = evalFilterRes(aFreqName,'1997-01-14', '0,1', 'HER', [10000 75536],'linhor,1.6,1,linhor,0.6,0.3,ellihor,0.6,0.3');

%---------------------------------------------------------------------------
figure(1); 
aGeom.aPiLeft = 17.51; aGeom.aPiRight = 17.69; plotsPolarizRatioTilt('1996-12-07','HER', aGeom, 0,0.97,' ','HER, North','HER, East');
aGeom.aPiLeft = 15.78; aGeom.aPiRight = 15.94; plotsPolarizRatioTilt('1997-01-04','HER', aGeom, 1,0.65,' ','HER, North','HER, East');
aGeom.aPiLeft = 11.94; aGeom.aPiRight = 12.09; plotsPolarizRatioTilt('1997-01-14','HER', aGeom, 1,0.33,gwlGetNotation('TIME','hours'),'HER, North','HER, East');

%---------------------------------------------------------------------------
figure(2); 
aFilterName1 = 'Spectral bandpass filter';
aFilterName2 = 'Total horizontal filter';
plotsFilterRes(aTime,aSigBand1,aSigPolarLin1+aSigPolarElli1,1.5,0.91,' ','1996-12-07, HER, North','HER, East',1,17.6,aFilterName1,aFilterName2);
plotsFilterRes(aTime,aSigBand2,aSigPolarLin2+aSigPolarElli2,1.5,0.61,' ','1997-01-04, HER, North','HER, East',0,15.9,'','');
plotsFilterRes(aTime,aSigBand3,aSigPolarLin3+aSigPolarElli3,1.5,0.31,gwlGetNotation('TIME','hours'),'1997-01-14, HER, North','HER, East',0,12.1,'','');

%---------------------------------------------------------------------------
figure(3); 
aFilterName1 = 'Linear horizontal filter';
aFilterName2 = 'Elliptical horizontal filter';
plotsFilterRes(aTime,aSigPolarLin1,aSigPolarElli1,1.5,0.91,' ','1996-12-07, HER, North','HER, East',1,17.6,aFilterName1,aFilterName2);
plotsFilterRes(aTime,aSigPolarLin2,aSigPolarElli2,1.5,0.61,' ','1997-01-04, HER, North','HER, East',0,15.9,'','');
plotsFilterRes(aTime,aSigPolarLin3,aSigPolarElli3,1.5,0.31,gwlGetNotation('TIME','hours'),'1997-01-14, HER, North','HER, East',0,12.1,'','');

%---------------------------------------------------------------------------
pause(0.00001);
delete('*.dat');
clear all;
print -f1 -r600 -depsc filterHERFig1;
print -f2 -r600 -depsc filterHERFig2;
print -f3 -r600 -depsc filterHERFig3;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function [aTime, aSigBand, aSigPolarLin, aSigPolarElli] = evalFilterRes(aFreqName,aFileName, aChann, aChanNot, aTimeInt, aFilter)
aNotation = strcat(aFileName,aChanNot);
aSignalName = strcat(aNotation,'sig1.dat');
aSpectrName = strcat(aNotation,'cwt1.dat');
gwlSignalRead(2,strcat(aFileName,'.asc'),'func',['--format=ASCII --chan=' aChann ' --tmin=' num2str(aTimeInt(1)) ' --tmax=' num2str(aTimeInt(2)) ' --istime --resample=16'],aSignalName,aNotation);
gwlExec('gwlET2DFilter',[' --infile=' aSignalName ' --outfile=' strcat(aNotation,'fil.dat') ' --filter=' aFilter ' --name="filtered signal" --type=complex --wavelet=cauchy --wavpar=10 --freq=' aFreqName]);
fid = fopen(strcat(aNotation,'fil(1).dat'),'r'); [aTime,aSigBand]=gwlReadSignal(fid); fclose(fid);
fid = fopen(strcat(aNotation,'fil(2).dat'),'r'); [aTime,aSigPolarLin]=gwlReadSignal(fid); fclose(fid);
fid = fopen(strcat(aNotation,'fil(3).dat'),'r'); [aTime,aSigPolarElli]=gwlReadSignal(fid); fclose(fid);

%---------------------------------------------------------------------------
function plotsFilterRes(aTime,aSigBand,aSigPolar,aMax,aY,aXlab,aY1lab,aY2lab,aColBar,aPi2x,aFilterName1,aFilterName2)
aTime = aTime./3600;
aWight = 0.12;
aPi2y = 1.4;

aX = 0.05;
aTitle = '';
if(aColBar == 1)
    aTitle = aFilterName1;
end;
gwlPlotFunction(aTime,real(aSigBand),aX,aY-aWight,0.43,aWight,min(aTime),max(aTime),-aMax,aMax,'',' ',aY1lab,aTitle);
    gwlText(aPi2x, aPi2y, 'Pi 2');
gwlPlotFunction(aTime,imag(aSigBand),aX,aY-2*aWight-0.008,0.43,aWight,min(aTime),max(aTime),-aMax,aMax,aXlab,' ',aY2lab);
    gwlText(aPi2x, aPi2y, 'Pi 2');
    
aX = 0.50;
aTitle = '';
if(aColBar == 1)
    aTitle = aFilterName2;
end;
gwlPlotFunction(aTime,real(aSigPolar),aX,aY-aWight,0.43,aWight,min(aTime),max(aTime),-aMax,aMax,'',' ',aY1lab,aTitle);
    gwlText(aPi2x, aPi2y, 'Pi 2');
gwlPlotFunction(aTime,imag(aSigPolar),aX,aY-2*aWight-0.008,0.43,aWight,min(aTime),max(aTime),-aMax,aMax,aXlab,' ',aY2lab);
    gwlText(aPi2x, aPi2y, 'Pi 2');
    

%---------------------------------------------------------------------------
function plotsPolarizRatioTilt(aFileName,aChanNot,aGeom,aColBar,aY,aXlab,aY1lab,aY2lab)
aSignalName = strcat(aFileName,aChanNot,'sig.dat');
aElliparName = strcat(aFileName,aChanNot,'par.dat');
fid = fopen(aSignalName,'r');    [aTime,aSig] = gwlReadSignal(fid);    fclose(fid);  
fid = fopen(aElliparName,'r');   [aTime,aFreq,aCWT] = gwlReadSpectrum(fid);   fclose(fid);  
aTime = aTime./3600;
aInd1 = 1+aGeom.aCat;
aInd2 = length(aTime)-aGeom.aCat;
aRatio = aCWT(:,:,2);
aRatio(1,aInd1+1) = 0;
aRatio(1,aInd1+2) = 1;
aTilt = aCWT(:,:,3);
aTilt(1,aInd1+1) = 90;
aTilt(1,aInd1+2) = -90;
aPuls1 = 150;
aPuls2 = 250;
aX = 0.04;
if(aColBar == 0)
    aLength = 0.463;
else  
    aLength = 0.395;
end;  
gwlPlotFunction(aTime,real(aSig),aX,aY-0.1,0.395,0.08,aTime(aInd1),aTime(aInd2),min(real(aSig)),max(real(aSig)),'',' ',aY1lab,aFileName);
    line([aGeom.aPiLeft,aGeom.aPiLeft],[min(real(aSig)),max(real(aSig))],'Color',gwlGetColor(2));
    line([aGeom.aPiRight,aGeom.aPiRight],[min(real(aSig)),max(real(aSig))],'Color',gwlGetColor(2));
gwlPlotImage(aTime(aInd1:aInd2),aFreq,aRatio(:,aInd1:aInd2),aX,aY-0.268,aLength,0.165,aXlab,gwlGetNotation('FREQ'),gwlGetNotation('EPAR','RATIO','WG'));
    line([aGeom.aPiLeft,aGeom.aPiLeft],[min(aFreq),max(aFreq)],'Color',gwlGetColor(2));
    line([aGeom.aPiRight,aGeom.aPiRight],[min(aFreq),max(aFreq)],'Color',gwlGetColor(2));
    set(gca,'YAxisLocation','right');
    if(aColBar == 0)
        set(gca,'YTickLabel',{});
        ylabel('');
        clb = colorbar;
        set(clb,'FontSize',7);
    end;  
gwlPlotFunction(aTime,imag(aSig),aX+0.48,aY-0.1,0.395,0.08,aTime(aInd1),aTime(aInd2),min(imag(aSig)),max(imag(aSig)),'',' ',aY2lab,aFileName);
    line([aGeom.aPiLeft,aGeom.aPiLeft],[min(imag(aSig)),max(imag(aSig))],'Color',gwlGetColor(2));
    line([aGeom.aPiRight,aGeom.aPiRight],[min(imag(aSig)),max(imag(aSig))],'Color',gwlGetColor(2));
gwlPlotImage(aTime(aInd1:aInd2),aFreq,aTilt(:,aInd1:aInd2),aX+0.48,aY-0.268,aLength,0.165,aXlab,gwlGetNotation('FREQ'),gwlGetNotation('EPAR','TILT','WG'));
    line([aGeom.aPiLeft,aGeom.aPiLeft],[min(aFreq),max(aFreq)],'Color',gwlGetColor(2));
    line([aGeom.aPiRight,aGeom.aPiRight],[min(aFreq),max(aFreq)],'Color',gwlGetColor(2));
    set(gca,'YAxisLocation','right');
    if(aColBar == 0)
        set(gca,'YTickLabel',{});
        ylabel('');
        clb = colorbar;
        set(clb,'FontSize',7);
    end;  
    