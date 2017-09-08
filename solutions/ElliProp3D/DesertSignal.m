function DesertSignal()
% The real seismograms (Fig. 1) used in example are a three-component record from 
% an explosive-source experiment aimed at imaging the Dead Sea Transform in the 
% Middle East (DESERT Group 2000. Multinational geoscientific research kicks off 
% in the Middle East. EOS, 81(50), 609, 616-617.). The horizontal seismograms were 
% rotated into radial and transverse components using the shot/receiver azimuth 
% obtained from the survey geometry.
% 
% [1] M.S.Diallo, M.Kulesh, M.Holschneider and F.Scherbaum. Instantaneous 
%     polarization attributes in the time-frequency domain and wave field separation 
%     // Preprint Series DFG SPP 1114, University of Bremen. Preprint 57 (2004).
% [2] M.S.Diallo, M.Kulesh, M.Holschneider and F.Scherbaum. Instantaneous 
%     polarization attributes in the time-frequency domain and wave field separation 
%     // Geophysical Prospecting. V. 53. No. 5. P. 723-731 (2005).
% [3] Michail A. Kulesh, Matthias Holschneider, Mamadou S. Diallo, Kristina I. 
%     Kurennaya. Elliptic properties of elastic surface waves in wavelet domain // 
%     Proceedings of the International Summer School "Advanced Problems in Mechanics" 
%     (Russia, St. Petersburg, June 28-July 5, 2005). P. 361-366 (2005).
% [4] Diallo M.S., Kulesh M., Holschneider M. and Scherbaum F. Instantaneous 
%     polarization attributes in the time-frequency domain: application to wave field 
%     separation // Eos Trans. AGU, 85(47), Fall Meet. Suppl., Abstract S31B-1063 
%     (2004).
% [5] Diallo M.S., Kurennaya K., Kulesh M. and Holschneider M. Elliptic 
%     properties of surface elastic waves in wavelet domain (in Russian) // Book of 
%     abstracts of XIV Winter School on Continuous Media Mechanics (28 February - 3 
%     March 2005, Perm). P. 100.
% 
% FIGURE 1. Three-component real seismograms showing radial, transverse and 
% vertical components. The strong surface-wave arrivals can be observed between 3 
% s and 5 s.
% 
% FIGURE 2. Wavelet transform moduli of the different signals in Fig. 1. Strong 
% energy arrivals that can be associated with the Rayleigh-wave arrivals can be 
% observed on (a) and (c) between approximately 3.8 s and 5 s.
% 
% FIGURE 3. Instantaneous polarization parameters used for the wave-type 
% separation: (a) the reciprocal ellipticity; (b) the angle between the direction 
% of the planarity vector and the x-axis; (c) the angle between the direction of 
% the planarity vector and the y-axis; (d) the angle between the direction of the 
% planarity vector and the z-axis.
% 
% FIGURE 4. Example of wave-type separation applied to the three-component 
% seismograms in Fig. 1. (a) Application of the polarization filter to identify 
% the seismic events with pure polarization in (x - y). Note the low amplitudes on 
% the filtered vertical component. (b) Filtered threecomponent seismograms 
% resulting from the multiplication of the wavelet transforms of each component by 
% the corresponding normalized cosine direction of the planarity vector followed 
% by an inverse wavelet transform. This filtering has the effect of reducing the 
% influence of energy arrivals that come from directions not parallel to the 
% direction in which each component is most sensitive to particle motion.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aSignalName = 'signal.dat';
aFreqName = 'freq.dat';
aSpectrName = 'cwt.dat';
aElliparName = 'elli.dat';
aFilteredName = 'filtered.dat';
aTmin = 0.0;
aTmax = 10;
aYmax = 180;
aImax = 2002;
aTK = 0.01;

%---------------------------------------------------------------------------
[aTime,aSignal,aParSig] = gwlSignalRead(1,'DesertSignal.sta','seis','--format=STA --geo=1 --chan=1,2,0 --tmin=0 --tmax=9.8 --to2p --rot=0,1,-1.163',aSignalName,'Experimental signal');

gwlCreateAxis(128,0.01,10,'lin',aFreqName,'Frequency');
gwlCwt(1, aSignalName, aFreqName, 2, 'morlet', 1, aSpectrName);
[aTime,aFreq,aCwt] = gwlConvert('3','',aSpectrName);

gwlExec('gwlET3D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=morozov --filter=5 --tw=1']);
[aTimeE,aFreqE,aCwtE] = gwlConvert('4,15,16,17','',aElliparName);

gwlExec('gwlET3DFilter',[' --infile=' aSpectrName ' --outfile=' aFilteredName ' --elli=' aElliparName ' --filter=norm,ellixy,0.7,0']);
[aTimeInv,aSignalInv] = gwlIwt(1, aFilteredName, 'delta');

%---------------------------------------------------------------------------
figure(1);
aSignal = fliplr(aSignal)*aYmax;
gwlPlotSeis(aTime,aSignal,0.07,0.37,0.9,0.47,aTmin,aTmax,1,gwlGetNotation('TIME'),'',0,'');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1,gwlGetNotation('MSIG','T','v'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2,gwlGetNotation('MSIG','T','t'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3,gwlGetNotation('MSIG','T','r'));

%---------------------------------------------------------------------------
figure(2);
gwlPlotImage(aTime(1:aImax),aFreq,aCwt(:,1:aImax,1),0.1,0.64,0.8,0.240,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('MSIG','WABS','r')]);
gwlPlotImage(aTime(1:aImax),aFreq,aCwt(:,1:aImax,2),0.1,0.39,0.8,0.240,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('MSIG','WABS','t')]);
gwlPlotImage(aTime(1:aImax),aFreq,aCwt(:,1:aImax,3),0.1,0.14,0.8,0.240,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('MSIG','WABS','v')]);

%---------------------------------------------------------------------------
figure(3);
WgRatio = aCwtE(:,1:aImax,1);
WgTilt1 = aCwtE(:,1:aImax,2);
WgTilt2 = aCwtE(:,1:aImax,3);
WgTilt3 = aCwtE(:,1:aImax,4);
WgRatio(1,1)=1;
gwlPlotImage(aTimeE(1:aImax),aFreqE,WgRatio,0.07,0.76,0.9,0.22,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('EPAR','RATIO','WG')]);
    colorbar;
    set(gwlColorBar,'YTick',0:0.2:1);
    set(gwlColorBar,'YTickLabel',{0,0.2,0.4,0.6,0.8,1});
WgTilt1(1,1) = -pi/2;    
WgTilt1(1,2) = pi/2;    
gwlPlotImage(aTimeE(1:aImax),aFreqE,WgTilt1,0.07,0.53,0.9,0.22,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','TILT','WG','x')]);
    colorbar;
    set(gwlColorBar,'YTick',-1.5:0.75:1.5);
    set(gwlColorBar,'YTickLabel',{'-pi/2','-pi/4','0','pi/4','pi/2'});
WgTilt2(1,1) = -pi/2;    
WgTilt2(1,2) = pi/2;    
gwlPlotImage(aTimeE(1:aImax),aFreqE,WgTilt2,0.07,0.3,0.9,0.22,'',gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('EPAR','TILT','WG','y')]);
    colorbar;
    set(gwlColorBar,'YTick',-1.5:0.75:1.5);
    set(gwlColorBar,'YTickLabel',{'-pi/2','-pi/4','0','pi/4','pi/2'});
WgTilt3(1,1) = -pi/2;    
WgTilt3(1,2) = pi/2;    
gwlPlotImage(aTimeE(1:aImax),aFreqE,WgTilt3,0.07,0.07,0.9,0.22,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(d) ' gwlGetNotation('EPAR','TILT','WG','z')]);
    colorbar;
    set(gwlColorBar,'YTick',-1.5:0.75:1.5);
    set(gwlColorBar,'YTickLabel',{'-pi/2','-pi/4','0','pi/4','pi/2'});
colormap hsv;

%---------------------------------------------------------------------------
figure(4);
aInvR = fliplr(aSignalInv(:,4:6))*aYmax;
gwlPlotSeis(aTimeInv,aInvR,0.07,0.54,0.9,0.44,aTmin,aTmax,1,' ','',0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1,gwlGetNotation('MSIG','T','v','r'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2,gwlGetNotation('MSIG','T','t','r'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3,gwlGetNotation('MSIG','T','r','r'));
aInvTh = fliplr(aSignalInv(:,1:3))*aYmax;
gwlPlotSeis(aTimeInv,aInvTh,0.07,0.07,0.9,0.44,aTmin,aTmax,1,gwlGetNotation('TIME'),'',0,'(b)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1,gwlGetNotation('MSIG','T','v','\theta'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2,gwlGetNotation('MSIG','T','t','\theta'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3,gwlGetNotation('MSIG','T','r','\theta'));

%---------------------------------------------------------------------------
pause(0.00001);
delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aElliparName); delete(aFilteredName);
clear all;

print -f1 -r600 -depsc DesertSignalFig1;
print -f2 -r600 -depsc DesertSignalFig2;
print -f3 -r600 -depsc DesertSignalFig3;
print -f4 -r600 -depsc DesertSignalFig4;
  
