function ChaosSig() 
% To make some sensible comparisons to the 'conventional' nearest neighbours method, 
% in this example we perform some tests on a typical and complicated chaotic time 
% series such as a time series generated from the Duffing equation, from Mackey-Glass 
% delay-differential equation and chaotic Ikeda map. These time series do 
% not have any global amplitude changing or trend, therefore we can perform the 
% direct comparison between three metrics.
% 
% [1] M. Kulesh, et al., Adaptive metrics in the nearest neighbours 
%     method, Physica D (2007), doi:10.1016/j.physd.2007.08.019
% 
% FIGURE 1. (a) Chaotic time series based on the solution of Duffing equation.
% (b) Horizontal component of the time series, vertical line shows the prediction 
% start. (c) Zoom of predicted values for three different metrics
%     
% FIGURE 2. (a) Chaotic time series based on the solution of the Mackey-Glass delay 
% differential equation. (b) Horizontal component of the time series, vertical line 
% shows the prediction start. (c) Zoom of predicted values for three different metrics    
% 
% FIGURE 3. (a) Chaotic time series based on the Ikeda map. (b) Vertical component of 
% the time series, vertical line shows the prediction start. (c) Zoom of predicted values 
% for three different metrics

%----------------------------------------------------------------------------
path(path, '../../mshell');
% The Mackey-Glass DDE
% We get for analyse the X-component of Mackey-Glass equation.
% The embedding dimension suggested by Taken's Embedding Theorem for the Mackey-Glass series 
% is 5 (Kaplan, D., Glass, L., Understanding Nonlinear Dynamics, Springer-Verlag, NY, 1995).
% localCreateMackeyGlass('ChaosSigMackeyGlass2D.asc');

% We get for analyse the (Y+4) component of Ikeda equation
% localIkeda('ChaosSigIkeda2D.asc');

%---------------------------------------------------------------------------
figure(1); 
aPar.aNameSource = 'ChaosSigDuffing.asc';
aPar.aName2D = 'ChaosSigDuffing2D.asc';
aPar.aNameL1 = 'ChaosSigDuffingL1.dat';
aPar.aNameL2 = 'ChaosSigDuffingL2.dat';
aPar.aNameL3 = 'ChaosSigDuffingL3.dat';
aPar.a2Dmin = -3;
aPar.a2Dmax = 3;
aPar.aName2Dx = 'x(t)';
aPar.aName2Dy = 'y(t)';
aPar.aNameSerie = 'x(t)+x_0';
localPlotOneSerie(aPar);

%---------------------------------------------------------------------------
figure(2); 
aPar.aNameSource = 'ChaosSigMackeyGlass.asc';
aPar.aName2D = 'ChaosSigMackeyGlass2D.asc';
aPar.aNameL1 = 'ChaosSigMackeyGlassL1.dat';
aPar.aNameL2 = 'ChaosSigMackeyGlassL2.dat';
aPar.aNameL3 = 'ChaosSigMackeyGlassL3.dat';
aPar.a2Dmin = 0.3;
aPar.a2Dmax = 1.4;
aPar.aName2Dx = 'x(t)';
aPar.aName2Dy = 'x(t-\tau_0)';
aPar.aNameSerie = 'x(t)';
localPlotOneSerie(aPar);

%---------------------------------------------------------------------------
figure(3); 
aPar.aNameSource = 'ChaosSigIkeda.asc';
aPar.aName2D = 'ChaosSigIkeda2D.asc';
aPar.aNameL1 = 'ChaosSigIkedaL1.dat';
aPar.aNameL2 = 'ChaosSigIkedaL2.dat';
aPar.aNameL3 = 'ChaosSigIkedaL3.dat';
aPar.a2Dmin = -12;
aPar.a2Dmax = 12;
aPar.aName2Dx = 'x(t)';
aPar.aName2Dy = 'y(t)';
aPar.aNameSerie = 'y(t)+y_0';
localPlotOneSerie(aPar);

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc ChaosSigFig1;
print -f2 -r600 -depsc ChaosSigFig2;
print -f3 -r600 -depsc ChaosSigFig3;


%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function localPlotOneSerie(aPar) 
[aTime,aSignal] = gwlSignalRead(2,aPar.aName2D,'func','--format=ASCII --istime');
aSource = load(aPar.aNameSource,'-ascii');
aRes1   = load(aPar.aNameL1,'-ascii'); 
aRes2   = load(aPar.aNameL2,'-ascii'); 
aRes3   = load(aPar.aNameL3,'-ascii'); 
aStart  = length(aSource)-length(aRes1)+1;
aEnd    = length(aSource);
aPar
calcMAPE(aSource, aRes1, aStart, aEnd)
calcMAPE(aSource, aRes2, aStart, aEnd)
calcMAPE(aSource, aRes3, aStart, aEnd)

gwlPlotFunction(real(aSignal),imag(aSignal),0.07,0.35,0.3,0.4,aPar.a2Dmin,aPar.a2Dmax,aPar.a2Dmin,aPar.a2Dmax,aPar.aName2Dx,aPar.aName2Dy,'(a)');
    grid off;

aMin = min(aSource);
aMax = max(aSource);
aAmpl=aMin:(aMax-aMin)/100:aMax;
gwlPlotFunction(1:aEnd,aSource,0.4,0.57,0.57,0.2,1,aEnd,aMin,aMax,' ',aPar.aNameSerie,'(b)');
    hold on; 
    plot(aStart:aEnd,aRes3,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);
    hold off;
    line([aStart,aStart],[aMin,aMax],'Color',gwlGetColor(1));
    grid off;
aMin = min(aSource(aStart:aEnd));
aMax = max(aSource(aStart:aEnd));
gwlPlotFunction(aStart:aEnd,aSource(aStart:aEnd),0.4,0.33,0.57,0.2,aStart,aEnd,aMin,aMax*1.1,'Time',aPar.aNameSerie,'(c)');
    hold on; 
    plot(aStart:aEnd,aRes3,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);
    plot(aStart:aEnd,aRes1,'Color',gwlGetColor(1),'LineStyle','-.','LineWidth',1);
    plot(aStart:aEnd,aRes2,'Color',gwlGetColor(1),'LineStyle',':','LineWidth',2);
    hold off;
    legend('Source data','Prediction using L_A','Prediction using L_1','Prediction using L_2');
    grid off;
    
function localCreateMackeyGlass(aSourceName)
sol = dde23(@ddes,17,0.5,[0, 4096]);
t = linspace(20,4096,2048);
x = deval(sol,t);
y = deval(sol,t - 17);
fid1=fopen(aSourceName,'w');
for i=1:length(t)
    fprintf(fid1,'%4.3e   %6.5e   %6.5e\n',t(i),x(i),y(i));
end
fclose(fid1);

function dydt = ddes(t,y,Z)
dydt = 0.2*Z/(1 + Z^10) - 0.1*y;

function localIkeda(aSourceName)
[t,x,y] = ikeda(2048,0,1,0,1);
fid1=fopen(aSourceName,'w');
for i=1:length(t)
    fprintf(fid1,'%4.3e   %6.5e   %6.5e\n',t(i),x(i),y(i));
end
fclose(fid1);


function [at,x,y]=ikeda(n,level,mu,x0,y0)
% Simulation of the Ikeda map.
%    x'=1+mu(xcos(t)-ysin(t)
%    y'=mu(xsin(t)+ycos(t))
% x and y are the simulated time series.
% n is the number of the simulated points.
% level is the noise standard deviation divided by the standard deviation of
%   the noise-free time series. We assume Gaussian noise with zero mean.
% mu is the parameter.
% x0 is the initial value for x.
% y0 is the initial value for y.
% Reference:
% Ikeda K (1979): Multiple-valued stationary state and its instability of the
% transmitted light by a ring cavity system. Optics Communications 30: 257

% Initialize
t=0.4-6/(1+x0^2+y0^2);
x(1,1)=1+mu*(x0*cos(t)-y0*sin(t));
y(1,1)=mu*(x0*sin(t)+y0*cos(t));
at(1,1)=1;
% Simulate
for i=2:n
    t=0.4-6/(1+x(i-1,1)^2+y(i-1,1)^2);
    x(i,1)=1+mu*(x(i-1,1)*cos(t)-y(i-1,1)*sin(t));
    y(i,1)=mu*(x(i-1,1)*sin(t)+y(i-1,1)*cos(t));
    at(i,1)=i;
end
% Add normal white noise
x=x+randn(n,1)*level*std(x);
y=y+randn(n,1)*level*std(y);
