function SynthSigB()
% This example is performed on a synthetic 3-C signal (not based on a physical 
% model) containing arrivals with different types of polarization as shown in 
% Figure 1. Wave-group A simulates a stationary ellipse; wave-group B simulates a 
% rotating ellipse, and wave-group C is that of a linearly polarized event. Wave-
% group D is a stationary ellipse in 3D space, and wave-group E corresponds to 
% that of a rotating ellipsoid.
% 
% [1] M.S.Diallo, M.Kulesh, M.Holschneider, K.Kurennaya and F.Scherbaum 
%     Instantaneous polarization attributes based on adaptive covariance method // 
%     Preprint Series DFG SPP 1114, University of Bremen. Preprint 87 (2005).
% [2] M.S.Diallo, M.Kulesh, M.Holschneider, K.Kurennaya, F.Scherbaum Instantaneous 
%     polarization attributes based on an adaptive approximate covariance method // 
%     Geophysics. V. 71. No. 5. P. V99-V104 (2006).
% [3] Diallo M.S., Kurennaya K., Kulesh M. and Holschneider M Elliptic properties 
%     of surface elastic waves in wavelet domain (in Russian) // Book of abstracts of 
%     XIV Winter School on Continuous Media Mechanics (28 February - 3 March 2005, 
%     Perm). P. 100.
% 
% FIGURE 1. (a) Three synthetic signals simulating a 3-C record with elliptically 
% polarized wave groups contained in the (x-z), (x-y), and (y-z) planes. Hodograph 
% plots show the ellipses in the different polarization plane for wave groups (b) 
% A-C, (c) for wave group D, and (d) for wave group E.
% 
% FIGURE 2. 3D instantaneous-polarization attributes in the time domain: (a) 
% instantaneous major, minor, and second minor axes; (b) ellipticity ratio; (c) 
% adaptive time window. From the ellipticity ratios, one can identify clearly the 
% linearly polarized wave group [Rmed(t)/Rmax(t) close to zero] from those with 
% elliptical or ellipsoidal polarization.
% 
% FIGURE 3. 3D instantaneous-polarization attributes in the time-frequency domain: 
% instantaneous minor (a), major (b), and second minor (c) axes.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aSignalName = 'signal.dat';
aFreqName = 'freq.dat';
aSpectrName = 'cwt.dat';
aElliparName = 'elli.dat';
aTmin = 0.0;
aTmax = 5;
aYmax = 10;
aImax = 2002;

%---------------------------------------------------------------------------
[aTime,aSignal] = localPrepareSignal(aSignalName);
[aTime,M1rmin,M1rmax,M1rmidd,M1ratio,M1ratio2,M1AdaptTime] = localCalcElliSignal(aSignalName,1);
[aTime,M2rmin,M2rmax,M2rmidd,M2ratio,M2ratio2,M2AdaptTime] = localCalcElliSignal(aSignalName,2);
[aTime,M3rmin,M3rmax,M3rmidd,M3ratio,M3ratio2,M3AdaptTime] = localCalcElliSignal(aSignalName,3);

gwlCreateAxis(128,0.1,60,'lin',aFreqName,'Frequency');
gwlCwt(1, aSignalName, aFreqName, 1, 'cauchy', 5, aSpectrName);

gwlExec('gwlET3D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=acovar --filter=10 --tw=1']);
[aTimeE,aFreqE,aCwtE] = gwlConvert('2,1,3','',aElliparName);

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aTime,aSignal(:,1),0.07,0.83,0.9,0.15,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','x'),'(a)');
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotFunction(aTime,aSignal(:,2),0.07,0.68,0.9,0.15,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','y'),'');
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotFunction(aTime,aSignal(:,3),0.07,0.53,0.9,0.15,aTmin,aTmax,-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('MSIG','T','z'),'');
    gwlText(0.4,-6,'Event A');
    gwlText(1.4,-6,'Event B');
    gwlText(2.4,-6,'Event C');
    gwlText(3.4,-6,'Event D');
    gwlText(4.4,-6,'Event E');
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
subplot('position',[0.1,0.07,0.25,0.35]);
    aInd1 = 1;
    aInd2 = 1200;
    plot3(aSignal(aInd1:aInd2,1),aSignal(aInd1:aInd2,2),aSignal(aInd1:aInd2,3),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    gwlTitle('(b) Events A-C');
    gwlLabel('X','X');
    gwlLabel('Y','Y');
    gwlLabel('Z','Z');
    axis([-aYmax,aYmax,-aYmax/2,aYmax/2,-aYmax/10,aYmax/10]); 
    grid on;
subplot('position',[0.4,0.07,0.25,0.35]);
    aInd1 = 1200;
    aInd2 = 1600;
    plot3(aSignal(aInd1:aInd2,1),aSignal(aInd1:aInd2,2),aSignal(aInd1:aInd2,3),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    gwlTitle('(c) Event D');
    gwlLabel('X','X');
    gwlLabel('Y','Y');
    gwlLabel('Z','Z');
    axis([-aYmax,aYmax,-aYmax/2,aYmax/2,-aYmax/2.5,aYmax/2.5]); 
    grid on;
subplot('position',[0.7,0.07,0.25,0.35]);
    aInd1 = 1600;
    aInd2 = 2000;
    plot3(aSignal(aInd1:aInd2,1),aSignal(aInd1:aInd2,2),aSignal(aInd1:aInd2,3),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    gwlTitle('(d) Event E');
    gwlLabel('X','X');
    gwlLabel('Y','Y');
    gwlLabel('Z','Z');
    axis([-aYmax,aYmax,-aYmax/2,aYmax/2,-aYmax/2.5,aYmax/2.5]); 
    grid on;
 
%---------------------------------------------------------------------------
figure(2);
gwlPlotFunction(    aTime,M1rmin,0.07,0.67,0.9,0.22,aTmin,aTmax,0,aYmax,'',[gwlGetNotation('EPAR','RMAX','T') ', ' gwlGetNotation('EPAR','RMED','T') ', ' gwlGetNotation('EPAR','RMIN','T')],'(a)');
    hold on;   plot(aTime,M2rmin,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M3rmin,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M1rmax,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M2rmax,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M3rmax,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M1rmidd,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M2rmidd,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M3rmidd,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    gwlText(0.6,8.5,gwlGetNotation('EPAR','RMAX','T'));
    gwlText(0.42,3.0,gwlGetNotation('EPAR','RMED','T'));
    gwlText(4.42,2.5,gwlGetNotation('EPAR','RMIN','T'));
gwlPlotFunction(    aTime,M1ratio,0.07,0.436,0.9,0.22,aTmin,aTmax,0,1,'',gwlGetNotation('EPAR','RATIO','T'),'(b)');
    hold on;   plot(aTime,M2ratio,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M3ratio,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M2ratio2,'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
    hold on;   plot(aTime,M1ratio2,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
gwlPlotFunction(    aTime,M1AdaptTime,0.07,0.203,0.9,0.22,aTmin,aTmax,0,1.5,gwlGetNotation('TIME'),gwlGetNotation('EPAR','ATIME','T'),'(c)');

%---------------------------------------------------------------------------
figure(3);
gwlPlotImage(aTimeE(1:aImax),aFreqE,aCwtE(:,1:aImax,3),0.07,0.70,0.9,0.27,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('EPAR','RMIN','WG')]);
    colorbar;
gwlPlotImage(aTimeE(1:aImax),aFreqE,aCwtE(:,1:aImax,2),0.07,0.40,0.9,0.27,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','RMAX','WG')]);
    colorbar;
gwlPlotImage(aTimeE(1:aImax),aFreqE,aCwtE(:,1:aImax,1),0.07,0.10,0.9,0.27,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('EPAR','RMED','WG')]);
    colorbar;

%---------------------------------------------------------------------------
pause(0.00001);
delete(aSignalName); delete(aFreqName); delete(aSpectrName); delete(aElliparName);
clear all;
print -f1 -r600 -depsc SynthSigBFig1;
print -f2 -r600 -depsc SynthSigBFig2;
print -f3 -r600 -depsc SynthSigBFig3;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function [aTime,aSignal] = localPrepareSignal(aSourceName)
Size=2048;
dt=5/2048;
f0=10;
frot=f0/10;
gamma=14;
ex=(2*pi*f0/gamma)*(2*pi*f0/gamma);
for i=1:Size
    xt(i)=dt*i;
    % Stationary ellipse in 2D
    t1=dt*(i-0.1*Size);
    x1(i)=8*cos(2.0*pi*f0*xt(i))*(exp(-ex*t1*t1));
    y1(i)=4*sin(2.0*pi*f0*xt(i))*(exp(-ex*t1*t1));
    z1(i)=0;
    % Rotating ellipse in 2D space
    t3=dt*(i-0.3*Size);
    x3s(i)=8*cos(2.0*pi*f0*xt(i))*(exp(-ex*t3*t3));
    y3s(i)=4*sin(2.0*pi*f0*xt(i))*(exp(-ex*t3*t3));
    z3s(i)=0;
    x3(i)=+cos(2*pi*frot*xt(i))*x3s(i)+sin(2*pi*frot*xt(i))*y3s(i);
    y3(i)=-sin(2*pi*frot*xt(i))*x3s(i)+cos(2*pi*frot*xt(i))*y3s(i);
    z3(i)=0;
    % linear polarized signal in 2D space 
    t4=dt*(i-0.5*Size);
    x4(i)=6*cos(3.0*pi*f0*xt(i))*(exp(-ex*t4*t4));
    y4(i)=0.1*sin(3.0*pi*f0*xt(i))*(exp(-ex*t4*t4));
    z4(i)=0;
    % stationary ellipse in 3D space
    t2=dt*(i-0.7*Size);
    x2(i)=8*cos(3.0*pi*f0*xt(i))*(exp(-ex*t2*t2));
    y2(i)=4*sin(3.0*pi*f0*xt(i))*(exp(-ex*t2*t2));
    z2(i)=0.3*x2(i)+0.8*y2(i);
    % Ellipsoid  in 3D space 
    t5=dt*(i-0.9*Size);
    x5s(i)=8*cos(4.0*pi*f0*xt(i))*(exp(-ex*t5*t5));
    y5s(i)=4*sin(4.0*pi*f0*xt(i))*(exp(-ex*t5*t5));
    x5(i)=+cos(2*pi*frot*xt(i))*x5s(i)+sin(2*pi*frot*xt(i))*y5s(i);
    y5(i)=-sin(2*pi*frot*xt(i))*x5s(i)+cos(2*pi*frot*xt(i))*y5s(i);
    z5(i)=0.3*x5s(i)+0.8*y5s(i);
end;
x=x1+x2+x3+x4+x5;
y=y1+y2+y3+y4+y5;
z=z1+z2+z3+z4+z5;
fid1=fopen(aSourceName,'w');
for i=1:Size
    fprintf(fid1,'%f %+4.3e %+4.3e %+4.3e\n',xt(i),x(i),y(i),z(i));
end
fclose(fid1);
[aTime,aSignal] = gwlSignalRead(1,aSourceName,'seis','--format=ASCII --istime',aSourceName,'Synthetic signal');

%---------------------------------------------------------------------------
function [aTime,aRmin,aRmax,aRmidd,aRatio,aRatios,aAtime] = localCalcElliSignal(aSourceName,aType)
aTmpName1 = tempname;
if(aType == 1)
    gwlExec('gwlET3D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=acovar --tw=2']);
end;
if(aType == 2)
    gwlExec('gwlET3D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=scovar --tw=44']);
end;
if(aType == 3)
    gwlExec('gwlET3D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=morozov']);
end;
gwlExec('gwlConvert',[' --infile=' aTmpName1 ' --outfile=' aTmpName1 ' --outtype=1 --comp=2,1,3,4,6,7,8,9 --nomess']);
[aTime,aRmin,aRmax,aRmidd,aRatio,aRatios,aWx,aWy,aWz] = textread(aTmpName1,'%f %f %f %f %f %f %f %f %f');
delete(aTmpName1);
if(aType == 1)
    aAtime = 6*pi./(aWx+aWy+aWz);    
else
    aAtime(1:length(aTime)) = 0;    
end;
