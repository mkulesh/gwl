function NorthernBavaria()
% We illustrate the importance of the application of the method developed in this 
% contribution to three-component broad-band data. Here we use a set of high 
% quality seismograms recorded at the GRA1 station of the GRF array in northern 
% Bavaria (Germany). The event [origin time: 13-mar-1989 13:02, 
% lat: 50.7198N, lon: 9.9112E, h = 0.7 km, ML = 5.6, On-line Bulletin (2001)] was 
% triggered by the collapse of a salt mine in Thueringen in 1989 in a distance of 
% 148 km and a backazimuth of 321.2. against North from station GRA1. Fig. 1 shows 
% also the wavelet transforms of the different components, respectively. The 
% seismograms start with the P-wave onset (about 165 s relative time, time axis 
% starts at 1989 March 13, 13:00), which at this distance is composed of the head 
% wave Pn and the direct P wave Pg in regional velocity models used for 
% localization.
%
% [1] M. Kulesh, M. S. Diallo, M. Holschneider, K. Kurennaya, F. Krueger, M. 
%     Ohrnberger, F. Scherbaum. Polarization analysis in wavelet domain based on the 
%     adaptive covariance method // Preprint Series DFG SPP 1114, University of 
%     Bremen. Preprint 137 (2005).
% [2] Mamadou S. Diallo, Michail Kulesh, Matthias Holschneider, Kristina Kurrenaya 
%     and Frank Scherbaum. Estimating polarization attributes with an adaptive 
%     covariance method in the wavelet domain // SEG Technical Program Expanded 
%     Abstracts. V. 24. P. 1014-1017 (2005).
% [3] M. Kulesh, M. S. Diallo, M. Holschneider, K. Kurennaya, F. Krueger, M. Ohrnberger 
%     and F. Scherbaum. Polarization analysis in the wavelet domain based on the adaptive 
%     covariance method // Geophys. J. Int, Vol. 170, N. 2, P. 667-679 (2007). 
%     doi: 10.1111/j.1365-246X.2007.03417.x
% 
% FIGURE 1. Three-component seismograms recorded at the GRA1 stations, generated 
% by a salt mine collapse in Thuringen, Germany in 1989 and the continuous wavelet 
% transforms of this seismograms performed componentwise: (a) east component, (b) 
% north component and (c) vertical component.
% 
% FIGURE 2. A filtration based on the signed ellipticity ratio: (a) input data, 
% (b) positive part of the signed ellipticiy ratio in time–frequency plot and (c) 
% extracted prograde mutlicomponent signal with the filtration condition 
% WgRs(t,f)>0.15.
% 
% FIGURE 3. A filtration based on the signed ellipticity ratio: (a) input data, 
% (b) negative part of the signed ellipticiy ratio in time–frequency plot and (c) 
% extracted retrograde mutlicomponent signal with the filtration condition 
% WgRs(t,f)<-0.15.
% 
% FIGURE 4. A filtration based on the azimuth: (a) input data, (b) deviation of 
% the particle motions from the theoretically expected direction A0 = 128.8. and 
% (c) extracted mutlicomponent signal with the filtration conditions (WgA(t,f)-A0) 
% \in [-20, 20] and [70, -70] for the red line and (WgA(t,f)-A0) \in [25, 65] and 
% [-65, -25] for the black line.
% 
% FIGURE 5. A filtration based on the dip attribute: (a) input data, (b) the dip 
% attribute given in time-frequency plot and (c) polarization filtered results of 
% the input data to extract part of the multicomponent signal arriving at an 
% apparent dip angle range defined by the condition WgD(t,f) \in [70, 90].

%---------------------------------------------------------------------------
path(path, '../../mshell');

aSignalName = 'NorthernBavariaSig.dat';
aFreqName = 'freq.dat';
aSpectrName = 'NorthernBavariaCwt.dat';
aElliparName = 'NorthernBavariaElli.dat';
aFilteredName = 'filtered.dat';
[aTime,aSignal] = gwlSignalRead(1,'NorthernBavaria.asc','seis','--istime --chan=2,1,0 --tmin=150 --tmax=250 --to2p',aSignalName,'3C experimental seismogram');
gwlCreateAxis(128,0.0001,1,'lin',aFreqName,'Frequency');
gwlCwt(1, aSignalName, aFreqName, 2, 'morlet', 1.5, aSpectrName);
[aTime,aFreq,aCwt] = gwlConvert('3','',aSpectrName);
gwlExec('gwlET3D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=acovar --filter=1 --tw=3']);
[aTimeE,aFreqE,aCwtE] = gwlConvert('4,18,20,19','',aElliparName);
gwlExec('gwlET3DFilter',[' --infile=' aSpectrName ' --outfile=' aSpectrName ' --elli=' aElliparName ' --filter=prog,0.15,1,retro,0.15,1,azimuth,1.9,2.6,azimuth,-1.24,-0.54,azimuth,0.33,1.03,azimuth,-2.81,-2.11,dip,1.22,1.57']);
[aTimeInv,aSignalInv] = gwlIwt(1, aSpectrName, 'delta');
delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aElliparName);

%---------------------------------------------------------------------------
aYmax = 1/4000000;
aTmin = min(aTime);
aTmax = max(aTime);
aFmin = min(aFreq);
aFmax = max(aFreq);
aTK = 0.015;
aBand = 0.71;
aLogScale = 2000;
aCwt(:,:,:) = log10(aCwt(:,:,:)+aLogScale);
amax = max(max(max(aCwt)));
SFilter = sqrt(aCwt(:,:,1).^2+aCwt(:,:,2).^2+aCwt(:,:,3).^2);
WgBand =  SFilter./max(max(SFilter));
[aSizeF,aSizeT]=size(aCwt(:,:,1));
for i=1:aSizeF
    for j=1:aSizeT
        if WgBand(i,j)<aBand
            WgBand(i,j) = 0;    
        end;
    end;
end;
aSignal = fliplr(aSignal)*aYmax;
aPrograde = fliplr(aSignalInv(:,1:3)*aYmax);
aRegrade = fliplr(aSignalInv(:,4:6)*aYmax);
aAzimith1 = fliplr((aSignalInv(:,7:9)+aSignalInv(:,10:12))*aYmax);
aAzimith2 = fliplr((aSignalInv(:,13:15)+aSignalInv(:,16:18))*aYmax);
aDip = fliplr(aSignalInv(:,19:21)*aYmax);
clear SFilter;

%---------------------------------------------------------------------------
figure(1);
aHight1 = 0.2;
aHight2 = 0.1;
gwlPlotImage(aTime,aFreq,aCwt(:,:,1),0.07,0.09+2*(aHight1+aHight2),0.9,aHight1,'',gwlGetNotation('FREQ'),gwlGetNotation('MSIG','WABS','e'));
gwlPlotFunction(aTime,aSignal(:,3).*1.5,0.07,0.09+3*aHight1+2*aHight2,0.9,aHight2,aTmin,aTmax,-1,1,'','',['(a) ' gwlGetNotation('MSIG','T','e')]);
gwlPlotImage(aTime,aFreq,aCwt(:,:,2),0.07,0.08+1*(aHight1+aHight2),0.9,aHight1,'',gwlGetNotation('FREQ'),gwlGetNotation('MSIG','WABS','n'));
gwlPlotFunction(aTime,aSignal(:,2).*1.5,0.07,0.08+2*aHight1+aHight2,0.9,aHight2,aTmin,aTmax,-1,1,'','',['(b) ' gwlGetNotation('MSIG','T','n')]);
gwlPlotImage(aTime,aFreq,aCwt(:,:,3),0.07,0.07,0.9,aHight1,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),gwlGetNotation('MSIG','WABS','z'));
gwlPlotFunction(aTime,aSignal(:,1).*1.5,0.07,0.07+aHight1,0.9,aHight2,aTmin,aTmax,-1,1,'','',['(c) ' gwlGetNotation('MSIG','T','z')]);
colormap gwlgray;

%---------------------------------------------------------------------------
figure(2);
WgSignRatio = -aCwtE(:,:,2).*aCwtE(:,:,1);
for i=1:aSizeF
  for j=1:aSizeT
      if WgBand(i,j)==0
         WgSignRatio(i,j) = -0.05;    
      end;
      if WgSignRatio(i,j)<0
         WgSignRatio(i,j) = -0.05;    
      end;
  end;
end;
WgSignRatio(1,1)=1.05;
WgSignRatio(1,2)=-0.05;
gwlPlotSeis(aTime,aSignal,0.07,0.73,0.77,0.25,aTmin,aTmax,1,'','',0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.6,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.6,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.6,gwlGetNotation('MSIG','T','e'));
gwlPlotImage(aTime,aFreq,WgSignRatio,0.07,0.40,0.9,0.315,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','RATIOS','WG') '<0']);
    colorbar;
    set(gwlColorBar,'YTick',[0 0.2 0.4 0.6 0.8 1]);
    set(gwlColorBar,'YTickLabel',[0 -0.2 -0.4 -0.6 -0.8 -1]);
    set(gwlColorBar,'Box','off');
gwlPlotSeis(aTimeInv,aRegrade,0.07,0.07,0.77,0.31,aTmin,aTmax,1,gwlGetNotation('TIME'),'',0,['(c) ' gwlGetNotation('EPAR','RATIOS','WG') '<-0.15']);
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.5,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.5,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.5,gwlGetNotation('MSIG','T','e'));
colormap gwlgray;
clear WgSignRatio;

%---------------------------------------------------------------------------
figure(3);
WgSignRatio = aCwtE(:,:,2).*aCwtE(:,:,1);
for i=1:aSizeF
  for j=1:aSizeT
      if WgBand(i,j)==0
         WgSignRatio(i,j) = -0.05;    
      end;
      if WgSignRatio(i,j)<0
         WgSignRatio(i,j) = -0.05;    
      end;
  end;
end;
WgSignRatio(1,1)=1.05;
WgSignRatio(1,2)=-0.05;
gwlPlotSeis(aTime,aSignal,0.07,0.73,0.77,0.25,aTmin,aTmax,1,'','',0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.6,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.6,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.6,gwlGetNotation('MSIG','T','e'));
gwlPlotImage(aTime,aFreq,WgSignRatio,0.07,0.40,0.9,0.315,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','RATIOS','WG') '>0']);
    colorbar;
    set(gwlColorBar,'YTick',[0 0.2 0.4 0.6 0.8 1]);
    set(gwlColorBar,'YTickLabel',[0 0.2 0.4 0.6 0.8 1]);
    set(gwlColorBar,'Box','off');
gwlPlotSeis(aTimeInv,aPrograde,0.07,0.07,0.77,0.31,aTmin,aTmax,1,gwlGetNotation('TIME'),'',0,['(c) ' gwlGetNotation('EPAR','RATIOS','WG') '>0.15']);
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.5,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.5,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.5,gwlGetNotation('MSIG','T','e'));
colormap gwlgray;
clear WgSignRatio;

%---------------------------------------------------------------------------
figure(4);
WgAzimuth = 180/pi*aCwtE(:,:,3);
aMult = 1.8;
aAz0 = 128.8;
for i=1:aSizeF
    for j=1:aSizeT
        WgAzimuth(i,j) = WgAzimuth(i,j)-aAz0;
        if WgAzimuth(i,j) < -180
            WgAzimuth(i,j) = WgAzimuth(i,j) + 360; 
        end;
        if WgAzimuth(i,j) > 90
            WgAzimuth(i,j) = 180 - WgAzimuth(i,j);    
        end;
        if WgAzimuth(i,j) < -90
            WgAzimuth(i,j) = 180 + WgAzimuth(i,j);    
        end;
        if WgBand(i,j)==0
            WgAzimuth(i,j) = -95;    
        end;
    end;
end;
WgAzimuth(1,1)=-95;
WgAzimuth(1,1)=95;
gwlPlotSeis(aTime,aSignal,0.07,0.73,0.77,0.25,aTmin,aTmax,1,'','',0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.6,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.6,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.6,gwlGetNotation('MSIG','T','e'));
gwlPlotImage(aTime,aFreq,WgAzimuth,0.07,0.40,0.9,0.315,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','AZ','WG') '-' gwlGetNotation('EPAR','BAZ')]);
    colorbar;
    set(gwlColorBar,'YTick',[-90 -60 -30 0 30 60 90]);
    set(gwlColorBar,'YTickLabel',[-90 -60 -30 0 30 60 90]);
    set(gwlColorBar,'Box','off');
gwlPlotSeis(aTimeInv,aAzimith1*aMult,0.07,0.07,0.77,0.31,aTmin,aTmax,1,gwlGetNotation('TIME'),'',0,'(c)');
    gwlPlotSeisAdd(aTimeInv,aAzimith2*aMult,1,0,1);
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.5,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.5,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.5,gwlGetNotation('MSIG','T','e'));
    objsline = findobj(gca,'Type','line');
    set(objsline,'LineStyle','-');
colormap gwlhsv;

%---------------------------------------------------------------------------
figure(5);
aMult = 1.5;
WgDip = 180/pi*abs(aCwtE(:,:,4));
for i=1:aSizeF
    for j=1:aSizeT
        if WgBand(i,j)==0
            WgDip(i,j) = 0;    
        end;
    end;
end;
WgDip(1,1)=0;
WgDip(1,1)=90;
gwlPlotSeis(aTime,aSignal,0.07,0.73,0.77,0.25,aTmin,aTmax,1,'','',0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.6,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.6,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.6,gwlGetNotation('MSIG','T','e'));
gwlPlotImage(aTime,aFreq,WgDip,0.07,0.40,0.9,0.315,'',gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('EPAR','DIP','WG')]);
    colorbar;
    set(gwlColorBar,'YTick',[0 15 30 45 60 75 90]);
    set(gwlColorBar,'YTickLabel',[0 15 30 45 60 75 90]);
    set(gwlColorBar,'Box','off');
gwlPlotSeis(aTimeInv,aDip*aMult,0.07,0.07,0.77,0.31,aTmin,aTmax,1,gwlGetNotation('TIME'),'',0,'(c)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.5,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.5,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.5,gwlGetNotation('MSIG','T','e'));
colormap gwlgray;

%---------------------------------------------------------------------------
% Additional pictures from the presentation [2]
figure(6);
gwlPlotSeis(aTime,aSignal,0.07,0.73,0.77,0.25,aTmin,aTmax,1,'','',0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.3,gwlGetNotation('MSIG','T','z'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.3,gwlGetNotation('MSIG','T','n'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.3,gwlGetNotation('MSIG','T','e'));
gwlPlotImage(aTime,aFreq,aCwt(:,:,1),0.07,0.51,0.9,0.21,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('MSIG','WABS','e')]);
    colorbar;
gwlPlotImage(aTime,aFreq,aCwt(:,:,2),0.07,0.29,0.77,0.21,'',gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('MSIG','WABS','n')]);
gwlPlotImage(aTime,aFreq,aCwt(:,:,3),0.07,0.07,0.77,0.21,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(d) ' gwlGetNotation('MSIG','WABS','z')]);
colormap hsv;
 
%---------------------------------------------------------------------------
figure(7);
gwlPlotImage(aTime,aFreq,aCwtE(:,:,1).*WgBand,0.07,0.76,0.9,0.22,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('EPAR','RATIO','WG')]);
    colorbar;
gwlPlotImage(aTime,aFreq,aCwtE(:,:,1).*aCwtE(:,:,2).*WgBand,0.07,0.53,0.9,0.22,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','RATIOS','WG')]);
    colorbar;
gwlPlotImage(aTime,aFreq,WgAzimuth,0.07,0.30,0.9,0.22,'',gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('EPAR','DIP','WG')]);
    colorbar;
gwlPlotImage(aTime,aFreq,WgDip,0.07,0.07,0.9,0.22,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(d) ' gwlGetNotation('EPAR','AZ','WG')]);
    colorbar;
colormap hsv;

%---------------------------------------------------------------------------
pause(0.00001);
clear all;
  
print -f1 -r600 -depsc NorthernBavariaFig1;
print -f2 -r600 -depsc NorthernBavariaFig2;
print -f3 -r600 -depsc NorthernBavariaFig3;
print -f4 -r600 -depsc NorthernBavariaFig4;
print -f5 -r600 -depsc NorthernBavariaFig5;
 
