function ChileBolivia()
% In this example we tested the adaptive-covariance method (ACM) on the real 3-C 
% seismic record shown in Figure 1 and compared the result with the SCM. We 
% compared also the results of the covariancebased filter (Kanasewich, E. R., 
% 1981, Time sequence analysis in geophysics: University of Alberta Press.) on 
% this seismogram.
% 
% [1] M.S.Diallo, M.Kulesh, M.Holschneider, K.Kurennaya, F.Scherbaum Instantaneous 
%     polarization attributes based on an adaptive approximate covariance method // 
%     Geophysics. V. 71. No. 5. P. V99-V104 (2006).
% 
% [2] M. Kulesh, M. S. Diallo, M. Holschneider, K. Kurennaya, F. Krueger, M. 
%     Ohrnberger, F. Scherbaum. Polarization analysis in wavelet domain based on the 
%     adaptive covariance method // Preprint Series DFG SPP 1114, University of 
%     Bremen. Preprint 137 (2005).
% 
% FIGURE 1. 3-C seismogram originating from the June 13, 2005 earthquake in the 
% Chile- Bolivia border region recorded at the German Regional Seismic Network 
% (GRSN) station. The magnitude was determined at MW = 7.8.
% 
% FIGURE 2. Comparison of the polarization attributes obtained from the ACM with 
% those computed using SCM. (a) Major polarization axis, (b) second major 
% polarization axis, (c) minor polarization axis. While the twomethods agree well 
% on the average trend of curves, our method showsmore sensitivity to local 
% variations of the attributes.
% 
% FIGURE 3. Comparison of filtering results with the covariance-based filter. (a) 
% East component of the raw 3-C seismograms in Figure 1 with its filtered version 
% using the SCM and the ACM; (b) and (c) show the north and vertical components. 
% For the results from SCM, the window length T is given in several samples. For 
% the seismograms obtained with the ACM, the window length is indirectly 
% determined by the value N.
% 
% FIGURE 4. (a) Three-component seismograms for the MW = 7.8 earthquake of June 
% 13th, 2005 at the Chilean-Bolivian border recorded at the GRA1 station with 
% focal depth of 114 km.
% 
% FIGURE 5. The wavelet transforms of these seismograms performed componentwise. 
% We will focus our attention on the time windows for the Pdiff and SKSac phases.
% 
% FIGURE 6. Polarization analysis restricted to the time window corresponding to 
% the Pdiff phase arrival: (a) three-component seismogram for the considered time 
% window, (b) the wavelet transform of the vertical component, (c) the reciprocal 
% ellipticity in time-frequency domain. As expected, the reciprocal ellipticity 
% for this phase which is linearly polarized is very small. Interestingly, this 
% plot offers the possibility of directly estimating the frequency dependence of 
% the polarization.
% 
% FIGURE 7. Polarization analysis restricted to the time window corresponding to 
% the SKSac phase arrival: (a) three-component seismogram for the considered time 
% window, (b) the wavelet transform of the vertical component, (c) the reciprocal 
% ellipticity as function time and frequency estimated using R(t,f) and r(t,f), 
% (d) the minor reciprocal ellipticity estimated with help of rs(t,f) and r(t,f).

%---------------------------------------------------------------------------
path(path, '../../mshell');
aSignalName = 'signal.dat';
aFreqName = 'freq.dat';
aSpectrName = 'cwt.dat';
aElliparName = 'elli.dat';
gwlCreateAxis(128,0.00001,0.2,'lin',aFreqName,'Frequency');

%---------------------------------------------------------------------------
figure(1);
aTmin = 400;
aTmax = 2000;
aYmax = 1/60000;
[aTime,aSignal] = gwlSignalRead(1,'ChileBolivia.asc','seis','--format=ASCII --istime --resample=8 --chan=0,2,1',aSignalName,'3D experimental seismogram');
gwlPlotFunction(aTime,aSignal(:,1)*aYmax,0.07,0.63,0.9,0.15,aTmin,aTmax,-1,1,'',gwlGetNotation('MSIG','T','r'),'');
gwlPlotFunction(aTime,aSignal(:,2)*aYmax,0.07,0.48,0.9,0.15,aTmin,aTmax,-1,1,'',gwlGetNotation('MSIG','T','t'),'');
gwlPlotFunction(aTime,aSignal(:,3)*aYmax,0.07,0.33,0.9,0.15,aTmin,aTmax,-1,1,gwlGetNotation('TIME'),gwlGetNotation('MSIG','T','v'),'');

%---------------------------------------------------------------------------
figure(2);
aSTname = 'SCM';
aADname = 'ACM';
[aTime,M1rmin,M1rmax,M1rmidd,M1inv]=localCalcElliSignal(aSignalName,1);
[aTime,M2rmin,M2rmax,M2rmidd,M2inv]=localCalcElliSignal(aSignalName,2);
[aTime,M3rmin,M3rmax,M3rmidd,M3inv]=localCalcElliSignal(aSignalName,3);
[aTime,M4rmin,M4rmax,M4rmidd,M4inv]=localCalcElliSignal(aSignalName,4);
gwlPlotFunction(    aTime,M2rmax*aYmax,0.07,0.67,0.9,0.22,aTmin,aTmax,0,1,'',gwlGetNotation('EPAR','RMAX','T'),'(a)');
    hold on;   plot(aTime,M4rmax*aYmax,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
    legend(aADname,aSTname);
gwlPlotFunction(    aTime,M2rmin*aYmax,0.07,0.436,0.9,0.22,aTmin,aTmax,0,1,'',gwlGetNotation('EPAR','RMED','T'),'(b)');
    hold on;   plot(aTime,M4rmin*aYmax,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
gwlPlotFunction(    aTime,M2rmidd*aYmax,0.07,0.203,0.9,0.22,aTmin,aTmax,0,1,gwlGetNotation('TIME'),gwlGetNotation('EPAR','RMIN','T'),'(c)');
    hold on;   plot(aTime,M4rmidd*aYmax,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;

%---------------------------------------------------------------------------
figure(3);
aYmax = 1/80000;
Filterx(:,1)=M1inv(:,1); 
Filterx(:,2)=M2inv(:,1);
Filterx(:,3)=M3inv(:,1);
Filterx(:,4)=M4inv(:,1);
Filterx(:,5)=aSignal(:,1);
gwlPlotSeis(aTime,Filterx*aYmax,0.07,0.67,0.9,0.29,aTmin,aTmax,1,'',gwlGetNotation('MSIG','T','t'),0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*4.7,'Source');
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*3.7,strcat(aSTname,', n=50'));
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*2.7,strcat(aSTname,', n=10'));
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*1.7,strcat(aADname,', n=10'));
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*0.7,strcat(aADname,', n=2'));
Filtery(:,1)=M1inv(:,2); 
Filtery(:,2)=M2inv(:,2);
Filtery(:,3)=M3inv(:,2);
Filtery(:,4)=M4inv(:,2);
Filtery(:,5)=aSignal(:,2);
gwlPlotSeis(aTime,Filtery*aYmax,0.07,0.37,0.9,0.29,aTmin,aTmax,1,'',gwlGetNotation('MSIG','T','r'),0,'(b)');
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*4.7,'Source');
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*3.7,strcat(aSTname,', n=50'));
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*2.7,strcat(aSTname,', n=10'));
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*1.7,strcat(aADname,', n=10'));
    gwlText(aTmin+(aTmax-aTmin)*0.015,1*0.7,strcat(aADname,', n=2'));
Filterz(:,1)=M1inv(:,3); 
Filterz(:,2)=M2inv(:,3);
Filterz(:,3)=M3inv(:,3);
Filterz(:,4)=M4inv(:,3);
Filterz(:,5)=aSignal(:,3);
gwlPlotSeis(aTime,Filterz*aYmax,0.07,0.07,0.9,0.29,aTmin,aTmax,1,gwlGetNotation('TIME'),gwlGetNotation('MSIG','T','v'),0,'(c)');
    gwlText(aTmin+(aTmax-aTmin)*0.1,1*4.7,'Source');
    gwlText(aTmin+(aTmax-aTmin)*0.1,1*3.7,strcat(aSTname,', n=50'));
    gwlText(aTmin+(aTmax-aTmin)*0.1,1*2.7,strcat(aSTname,', n=10'));
    gwlText(aTmin+(aTmax-aTmin)*0.1,1*1.7,strcat(aADname,', n=10'));
    gwlText(aTmin+(aTmax-aTmin)*0.1,1*0.7,strcat(aADname,', n=2'));

%---------------------------------------------------------------------------
figure(4);
[aTimeWT,aSignalWT] = gwlSignalRead(1,'ChileBolivia.asc','seis','--format=ASCII --istime --resample=32 --chan=0,2,1',aSignalName,'3D experimental seismogram');
gwlCwt(1, aSignalName, aFreqName, 1, 'morlet', 1.5, aSpectrName);
[aTimeWT,aFreqWT,aCwt] = gwlConvert('3','',aSpectrName);
aTmin = min(aTimeWT);
aTmax = max(aTimeWT);
aYmax = 1/100000;
aTK = 0.012;
aSignalWT = fliplr(aSignalWT)*aYmax;
gwlPlotSeis(aTimeWT,aSignalWT,0.07,0.37,0.9,0.47,aTmin,aTmax,1,gwlGetNotation('TIME'),'',0,'');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1,gwlGetNotation('MSIG','T','v'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2,gwlGetNotation('MSIG','T','t'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3,gwlGetNotation('MSIG','T','r'));
    X1 = 480;  line([X1,X1],[0.4,1.6],'Color',gwlGetColor(2));  gwlText(X1+10,1.6,'Pdiff');
    X1 = 1100; line([X1,X1],[2.4,3.6],'Color',gwlGetColor(2));  gwlText(X1-160,3.6,'SKS');
    X1 = 1235; line([X1,X1],[2.4,3.6],'Color',gwlGetColor(2));  gwlText(X1+10,3.6,'Sdiff');

%---------------------------------------------------------------------------
figure(5);
aFmin = min(aFreqWT);
aFmax = max(aFreqWT);
amax = max(max(max(aCwt))); aCwt(1,1,1)=amax; aCwt(1,1,2)=amax; aCwt(1,1,3)=amax; 
gwlPlotImage(aTimeWT,aFreqWT,aCwt(:,:,1),0.1,0.64,0.8,0.240,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('MSIG','WABS','r')]);
gwlPlotImage(aTimeWT,aFreqWT,aCwt(:,:,2),0.1,0.39,0.8,0.240,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('MSIG','WABS','t')]);
gwlPlotImage(aTimeWT,aFreqWT,aCwt(:,:,3),0.1,0.14,0.8,0.240,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('MSIG','WABS','v')]);

%---------------------------------------------------------------------------
figure(6);
[aTimeP1,aSignalP1] = gwlSignalRead(1,'ChileBolivia.asc','seis','--format=ASCII --istime --chan=0,2,1 --tmin=450 --tmax=600 --resample=3 --to2p',aSignalName,'3D experimental seismogram');
gwlCwt(1, aSignalName, aFreqName, 1, 'morlet', 1.5, aSpectrName);
[aTimeP1,aFreqP1,aCwtP1] = gwlConvert('3','',aSpectrName);
gwlExec('gwlET3D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=acovar --filter=1 --tw=2']);
[aTimeP1,aFreqP1,aRatioP1] = gwlConvert('4','',aElliparName);
aTmin = min(aTimeP1);
aTmax = max(aTimeP1);
aFmin = min(aFreqP1);
aFmax = max(aFreqP1);
aYmax = 1/20000;
aSignalP1 = fliplr(aSignalP1)*aYmax;
gwlPlotSeis(aTimeP1,aSignalP1,0.07,0.73,0.77,0.25,aTmin,aTmax,1,'','',0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1,gwlGetNotation('MSIG','T','v'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2,gwlGetNotation('MSIG','T','t'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3,gwlGetNotation('MSIG','T','r'));
gwlPlotImage(aTimeP1,aFreqP1,aCwtP1(:,:,3),0.07,0.51,0.9,0.21,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('MSIG','WABS','v')]);
    colorbar;
gwlPlotImage(aTimeP1,aFreqP1,aRatioP1,0.07,0.29,0.9,0.21,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('EPAR','RATIO','WG')]);
    colorbar;

%---------------------------------------------------------------------------
figure(7);
[aTimeP2,aSignalP2] = gwlSignalRead(1,'ChileBolivia.asc','seis','--format=ASCII --istime --chan=0,2,1 --tmin=1100 --tmax=1151 --to2p',aSignalName,'3D experimental seismogram');
gwlCwt(1, aSignalName, aFreqName, 2, 'morlet', 1.5, aSpectrName);
[aTimeP2,aFreqP2,aCwtP2] = gwlConvert('3','',aSpectrName);
gwlExec('gwlET3D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=acovar --filter=1 --tw=2']);
[aTimeP2,aFreqP2,aRatioP2] = gwlConvert('4,5','',aElliparName);
aTmin = min(aTimeP2);
aTmax = max(aTimeP2);
aFmin = min(aFreqP2);
aFmax = max(aFreqP2);
aYmax = 1/100000;
aSignalP2 = fliplr(aSignalP2)*aYmax;
gwlPlotSeis(aTimeP2,aSignalP2,0.07,0.73,0.77,0.25,aTmin,aTmax,1,'','',0,'(a)');
    gwlText(aTmin+(aTmax-aTmin)*aTK,1,gwlGetNotation('MSIG','T','v'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2,gwlGetNotation('MSIG','T','t'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3,gwlGetNotation('MSIG','T','r'));
gwlPlotImage(aTimeP2,aFreqP2,aCwtP2(:,:,1),0.07,0.51,0.9,0.21,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('MSIG','WABS','r')]);
    colorbar;
gwlPlotImage(aTimeP2,aFreqP2,aRatioP2(:,:,1),0.07,0.29,0.9,0.21,'',gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('EPAR','RATIO','WG')]);
    colorbar;
gwlPlotImage(aTimeP2,aFreqP2,aRatioP2(:,:,2),0.07,0.07,0.9,0.21,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(d) ' gwlGetNotation('EPAR','RATIO1','WG')]);
    colorbar;

%---------------------------------------------------------------------------
pause(0.00001);
delete(aSignalName);  delete(aFreqName);  delete(aSpectrName);  delete(aElliparName);
clear all;

print -f1 -r600 -depsc ChileBoliviaFig1;
print -f2 -r600 -depsc ChileBoliviaFig2;
print -f3 -r600 -depsc ChileBoliviaFig3;
print -f4 -r600 -depsc ChileBoliviaFig4;
print -f5 -r600 -depsc ChileBoliviaFig5;
print -f6 -r600 -depsc ChileBoliviaFig6;
print -f7 -r600 -depsc ChileBoliviaFig7;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function [aTime,aRmin,aRmax,aRmidd,aSignalInv] = localCalcElliSignal(aSourceName,aType)
aTmpName1 = tempname;
aTmpName2 = tempname;
if(aType == 1)
    gwlExec('gwlET3D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=acovar --tw=2']);
end;
if(aType == 2)
    gwlExec('gwlET3D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=acovar --tw=10']);
end;
if(aType == 3)
    gwlExec('gwlET3D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=scovar --tw=10']);
end;
if(aType == 4)
    gwlExec('gwlET3D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=scovar --tw=50']);
end;
gwlExec('gwlET3DFilter',[' --infile=' aSourceName ' --outfile=' aTmpName2 ' --elli=' aTmpName1 ' --nomess']);
fid = fopen(aTmpName2,'r');     [aTime,aSignalInv]=gwlReadSignal(fid);  fclose(fid);
gwlExec('gwlConvert',[' --infile=' aTmpName1 ' --outfile=' aTmpName1 ' --outtype=1 --comp=2,1,3 --nomess']);
[aTime,aRmin,aRmax,aRmidd] = textread(aTmpName1,'%f %f %f %f');
delete(aTmpName1);
delete(aTmpName2);
