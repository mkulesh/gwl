function exmETacovar()
% exmETacovar(): In this example we show the effect of varying the parameter N
% in adaptive approximate covariance method.
% 
% [1] M.S.Diallo, M.Kulesh, M.Holschneider, K.Kurennaya, F.Scherbaum
%     Instantaneous polarization attributes based on an adaptive approximate 
%     covariance method // Geophysics. V. 71. No. 5. P. V99-V104 (2006).
% [2] M. Kulesh, M. S. Diallo, M. Holschneider, K. Kurennaya, F. Krueger, M. Ohrnberger, 
%     F. Scherbaum. Polarization analysis in wavelet domain based on the adaptive covariance method // 
%     Preprint Series DFG SPP 1114, University of Bremen. Preprint 137 (2005).
% 
% FIGURE 1. Illustration of the influence of increasing N on the accuracy of the approximation 
% for the signals. (a) Three synthetic signals Sx(t), Sy(t), Sz(t) (black solid lines) and the 
% corresponding approximations (dotted lines) for corresponding values of N = 1, 2, 3. (b, c, d) 
% Hodograph plots for each segment of the approximated signals for the different values of N. 
% For N = 1, the approximation is good up to 0.05 s, for N = 2 the approximation extends to 0.1 s 
% and for N = 3 the approximation covers the entirevtimewindow of the synthetic signals. Note how 
% the curves drift apart as N increases.
    
%---------------------------------------------------------------------------
path(path, '../../mshell');

aSignal = load('exmETacovar.asc','-ascii');
[N,N2]=size(aSignal);
aSamplPeriod = 2.442E-03;
aTime(1:N) = (0:N-1)*aSamplPeriod;
aYmax = -0.07;
aTmin = 0;
aTmax = 0.15;
aTK = 0.012;

Omz=120.633;   Hz=-6.323+4.888i;
Omy=111.714;   Hy=2.286+3.288i;
Omx=123.666;   Hx=-0.126-3.995i;

aApprox(:,2) = abs(Hx).*cos(Omx.*aTime(:) + angle(Hx));
aApprox(:,3) = abs(Hy).*cos(Omy.*aTime(:) + angle(Hy));
aApprox(:,1) = abs(Hz).*cos(Omz.*aTime(:) + angle(Hz));

%---------------------------------------------------------------------------
figure(1);
gwlPlotSeis(aTime,aSignal.*aYmax,0.07,0.53,0.9,0.45,aTmin,aTmax,1,gwlGetNotation('TIME'),'',0,'(a)');
    gwlPlotSeisAdd(aTime,aApprox.*aYmax,1,0,1);
    gwlText(aTmin+(aTmax-aTmin)*aTK,1.5,gwlGetNotation('MSIG','T','x'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,2.3,gwlGetNotation('MSIG','T','y'));
    gwlText(aTmin+(aTmax-aTmin)*aTK,3.3,gwlGetNotation('MSIG','T','z'));
subplot('position',[0.1,0.07,0.25,0.35]);
    N1=21;
    plot3(aApprox(1:N1,1),aApprox(1:N1,2),aApprox(1:N1,3),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    gwlTitle(['(b) n=1,   ' gwlGetNotation('EPAR','ATIME') '=' num2str(aTime(N1))]);
    gwlLabel('X','X');
    gwlLabel('Y','Y');
    gwlLabel('Z','Z');
    grid on;
    aTime(N1-1:N1+1)
    6*pi*1/(Omx+Omy+Omz)
subplot('position',[0.4,0.07,0.25,0.35]);
    N1=42;
    plot3(aApprox(1:N1,1),aApprox(1:N1,2),aApprox(1:N1,3),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    gwlTitle(['(c) n=2,   ' gwlGetNotation('EPAR','ATIME') '=' num2str(aTime(N1))]);
    gwlLabel('X','X');
    gwlLabel('Y','Y');
    gwlLabel('Z','Z');
    grid on;
    aTime(N1-1:N1+1)
    6*pi*2/(Omx+Omy+Omz)
subplot('position',[0.7,0.07,0.25,0.35]);
    N1=63;
    plot3(aApprox(1:N1,1),aApprox(1:N1,2),aApprox(1:N1,3),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    gwlTitle(['(d) n=3,   ' gwlGetNotation('EPAR','ATIME') '=' num2str(aTime(N1))]);
    gwlLabel('X','X');
    gwlLabel('Y','Y');
    gwlLabel('Z','Z');
    aTime(N1-1:N1+1)
    6*pi*3/(Omx+Omy+Omz)
    grid on;

%---------------------------------------------------------------------------
pause(0.00001);
clear all;

print -f1 -r600 -depsc exmETacovar;
