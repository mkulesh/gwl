function SynthSigA()
% In the following we show how to use the instantaneous attributes defined in the 
% wavelet domain for the purpose of separating different wave types. We simulate a 
% three-component record with seismic events consisting of the superposition of 
% elliptically and linearly polarized events. The simulation is not based on an 
% actual physical model but is rather considered as a simplified illustration of 
% the filtering principle. The actual simulated three-component record is shown in 
% Fig. 1. It consists of an elliptically polarized event in the (x - y) plane 
% between 1 s and 5 s, a second elliptically polarized event in the (x - z) plane 
% between 5 s and 9 s, a third elliptically polarized event in the (y - z) plane 
% between 9 s and 15 s and, finally, a mixture of linearly and elliptically 
% polarized events having different frequency content between 15 s and 20 s.
% 
% [1] M.S.Diallo, M.Kulesh, M.Holschneider and F.Scherbaum. Instantaneous 
%     polarization attributes in the time-frequency domain and wave field separation 
%     // Preprint Series DFG SPP 1114, University of Bremen. Preprint 57 (2004).
% 
% [2] M.S.Diallo, M.Kulesh, M.Holschneider and F.Scherbaum Instantaneous 
%     polarization attributes in the time-frequency domain and wave field separation 
%     // Geophysical Prospecting. V. 53. No. 5. P. 723-731 (2005).
% 
% FIGURE 1. (a) Three synthetic signals simulating a three-component record with 
% elliptically polarized events contained in the (x - z), (x - y) and (y - z) 
% planes. The late arrival on the three seismograms is a mixture of a linearly 
% polarized event and an elliptically polarized event with different frequency 
% content. (b) Hodogram plot showing the ellipses in the different polarization 
% planes. The dashed curves in Fig. 1(a) depict the filtering results, where from 
% the initial synthetic seismogram we have one distinguishable elliptical event 
% that extends from nearly 5 s to 10 s, and also the elliptical event in the time 
% window between 15 s and 20 s which is superimposed the linear event.
% 
% FIGURE 2. Instantaneous polarization attributes in the time-frequency domain: 
% (a) the ellipticity ratio for the different events in the multicomponent signal; 
% (b) the angle between the planarity vector and the x-axis; (c) the angle between 
% the planarity vector and the y-axis; (d) the angle between the planarity vector 
% and the z-axis. Note that for the first event, Theta_z(b, a) = 0; for the second 
% event, Theta_y(b, a) = 0; for the third event, Theta_x(b, a) = 0. Only the 
% linearly polarized event (high frequency) between 14 s and 20 s has non-zero 
% components for p in all directions. For the elliptically polarized event in that 
% time window Theta_y(b, a) = 0.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aSignalName = 'signal.dat';
aFreqName = 'freq.dat';
aSpectrName = 'cwt.dat';
aElliparName = 'elli.dat';
aFilteredName = 'filtered.dat';
aTmin = 0.0;
aTmax = 20;
aYmax = 2;
aImax = 2002;

%---------------------------------------------------------------------------
[aTime,aSignal] = localPrepareSignal(aSignalName);

gwlCreateAxis(128,0.01,5,'lin',aFreqName,'Frequency');
gwlCwt(1, aSignalName, aFreqName, 1, 'morlet', 1, aSpectrName);

gwlExec('gwlET3D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=morozov --filter=10 --tw=1']);
[aTimeE,aFreqE,aCwtE] = gwlConvert('4,15,16,17,19,20','--degree --modpi',aElliparName);

gwlExec('gwlET3DFilter',[' --infile=' aSpectrName ' --outfile=' aFilteredName ' --elli=' aElliparName ' --filter=ellixz,0.1,0.25']);
[aTimeInv,aSignalInv] = gwlIwt(1, aFilteredName, 'delta');

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aTime,aSignal(:,1),0.07,0.82,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','x'),'(a)');
    hold on;   plot(aTimeInv,aSignalInv(:,1),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;

gwlPlotFunction(aTime,aSignal(:,2),0.07,0.68,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','y'),'');
    hold on;   plot(aTimeInv,aSignalInv(:,2),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    
gwlPlotFunction(aTime,aSignal(:,3),0.07,0.54,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('MSIG','T','z'),'');
    hold on;   plot(aTimeInv,aSignalInv(:,3),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;

subplot('Position',[0.3,0.07,0.4,0.4]);
plot3(aSignal(:,1),aSignal(:,2),aSignal(:,3),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
gwlLabel('X','X');
gwlLabel('Y','Y');
gwlLabel('Z','Z');
axis([-aYmax,aYmax,-aYmax,aYmax,-aYmax,aYmax]); 
grid on;

%---------------------------------------------------------------------------
figure(2);
WgRatio = aCwtE(:,1:aImax,1);   WgRatio(1,1)=1;    
WgTilt1 = 180*aCwtE(:,1:aImax,2)/pi;   WgTilt1(1,1) = 0;    WgTilt1(1,2) = 95;    
WgTilt2 = 180*aCwtE(:,1:aImax,3)/pi;   WgTilt2(1,1) = 0;    WgTilt2(1,2) = 95;    
WgTilt3 = 180*aCwtE(:,1:aImax,4)/pi;   WgTilt3(1,1) = 0;    WgTilt3(1,2) = 95;    

gwlPlotImage(aTimeE(1:aImax),aFreqE,WgRatio,0.07,0.76,0.9,0.20,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('EPAR','RATIO','WG')]);
    colorbar;
    set(gwlColorBar,'YTick',0:0.2:1);
    set(gwlColorBar,'YTickLabel',{0,0.2,0.4,0.6,0.8,1});
gwlPlotImage(aTimeE(1:aImax),aFreqE,WgTilt1,0.07,0.53,0.9,0.22,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','TILT','WG','x')]);
    colorbar;
gwlPlotImage(aTimeE(1:aImax),aFreqE,WgTilt2,0.07,0.3,0.9,0.22,'',gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('EPAR','TILT','WG','y')]);
    colorbar;
gwlPlotImage(aTimeE(1:aImax),aFreqE,WgTilt3,0.07,0.07,0.9,0.22,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(d) ' gwlGetNotation('EPAR','TILT','WG','z')]);
    colorbar;

%---------------------------------------------------------------------------
figure(3);
[aSizeF,aSizeT]=size(aCwtE(:,:,6));
WgAzimuth = aCwtE(:,:,6);   WgAzimuth(1,1)=0;   WgAzimuth(1,2)=95;

gwlPlotImage(aTimeE,aFreqE,aCwtE(:,:,5),0.07,0.5,0.9,0.4,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('EPAR','DIP','WG')]);
    colorbar;
gwlPlotImage(aTimeE,aFreqE,abs(WgAzimuth),0.07,0.07,0.9,0.4,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','AZ','WG')]);
    colorbar;
colormap hsv;

%---------------------------------------------------------------------------
pause(0.00001);
delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aElliparName); delete(aFilteredName);
clear all;
print -f1 -r600 -depsc SynthSigAFig1;
print -f2 -r600 -depsc SynthSigAFig2;
print -f3 -r600 -depsc SynthSigAFig3;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function [aTime,aSignal] = localPrepareSignal(aSourceName)
dt=.01;
Nt=2*1024;
a=1.0;
for i=1:Nt
    t(i)=dt*(i-1);
    
    F1 = 2;    t1 = 3;
    x1(i)=cos(2*pi*F1*(t(i)-t1))*exp(-(t(i)-t1)*(t(i)-t1)/a);
    y1(i)=1.5*sin(2*pi*F1*(t(i)-t1))*exp(-(t(i)-t1)*(t(i)-t1)/a);
    z1(i)=0;
    
    F2 = 2;    t2 = 7.5;
    x2(i)=0.7*cos(2*pi*F2*(t(i)-t2))*exp(-(t(i)-t2)*(t(i)-t2)/a);
    y2(i)=0;
    z2(i)=1.3*sin(2*pi*F2*(t(i)-t2))*exp(-(t(i)-t2)*(t(i)-t2)/a);
    
    F3 = 2;    t3 = 12;
    x3(i)=0;
    y3(i)=0.5*cos(2*pi*F3*(t(i)-t3))*exp(-(t(i)-t3)*(t(i)-t3)/a);
    z3(i)=1.7*sin(2*pi*F3*(t(i)-t3))*exp(-(t(i)-t3)*(t(i)-t3)/a);
    
    F4 = 1;    t4 = 17;    F4lin = 3;    a=1.5;
    x4(i)=0.9*cos(2*pi*F4*(t(i)-t4))*exp(-(t(i)-t4)*(t(i)-t4)/a);
    y4(i)=0;
    z4(i)=1.5*sin(2*pi*F4*(t(i)-t4))*exp(-(t(i)-t4)*(t(i)-t4)/a);
    x4lin(i)=0.3*cos(2*pi*F4lin*(t(i)-t4))*exp(-(t(i)-t4)*(t(i)-t4)/a);
    y4lin(i)=0.6*cos(2*pi*F4lin*(t(i)-t4))*exp(-(t(i)-t4)*(t(i)-t4)/a);
    z4lin(i)=0.9*cos(2*pi*F4lin*(t(i)-t4))*exp(-(t(i)-t4)*(t(i)-t4)/a);
    
    x(i)=x1(i)+x2(i)+x3(i)+x4(i)+x4lin(i);
    y(i)=y1(i)+y2(i)+y3(i)+y4(i)+y4lin(i);
    z(i)=z1(i)+z2(i)+z3(i)+z4(i)+z4lin(i);
end;
fid1=fopen(aSourceName,'w');
for i=1:Nt
    fprintf(fid1,'%f %+4.3e %+4.3e %+4.3e\n',t(i),x(i),y(i),z(i));
end
fclose(fid1);
[aTime,aSignal] = gwlSignalRead(1,aSourceName,'seis','--format=ASCII --istime',aSourceName,'Synthetic signal');