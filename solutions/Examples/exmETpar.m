function exmETpar()
% exmETpar(): Schematic representation of an ellipse with its geometric parameters

%---------------------------------------------------------------------------
path(path, '../../mshell');

Z1mod = 2.7;
Z1arg = -0.2;
Z1 = Z1mod*exp(i*Z1arg);
Z2mod = 6.0;
Z2arg = 2.0;
Z2 = Z2mod*exp(i*Z2arg);
N = 1000;
for k=1:N
    aTime(k) = 2*pi*k/N;
    aElli(k) = Z1*exp(i*aTime(k))+Z2*exp(-i*aTime(k));
    aCirc(k) = 3*exp(i*aTime(k));
end;
aMax = 10;

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(real(aElli),imag(aElli),0.07,0.2,0.35,0.5,-aMax,aMax,-aMax,aMax,'','','');
    hold on;   plot(real(aElli),imag(aElli),'Color',gwlGetColor(0),'LineWidth',2);    hold off;
    line([0,0],[-aMax,aMax],'Color','black');
    line([-aMax,aMax],[0,0],'Color','black');
    aInd = 190; line([0,real(aElli(aInd))],[0,imag(aElli(aInd))],'Color','black');
    gwlText(real(2*aElli(aInd))/3,2*imag(aElli(aInd))/3-0.3,gwlGetNotation('EPAR','RMAX'));
    aInd = 930; line([0,real(aElli(aInd))],[0,imag(aElli(aInd))],'Color','black');
    gwlText(real(2*aElli(aInd))/3,2*imag(aElli(aInd))/3-0.6,gwlGetNotation('EPAR','RMIN'));
    hold on;   plot(real(aCirc(1:130)),imag(aCirc(1:130)),'-black','LineWidth',1);    hold off;
    aInd = 90;  gwlText(real(aCirc(aInd))+0.3,imag(aCirc(aInd))+0.3,gwlGetNotation('EPAR','TILT'));
    gwlText(9.0,-0.1,'x');
    gwlText(0.1,9.9,'z');
    set(gca,'Visible','Off');

gwlPlotFunction(aTime(1:N-50),real(aElli(1:N-50)),0.45,0.2,0.45,0.5,min(aTime),max(aTime),-aMax,aMax,'','','');
    hold on;   plot(aTime,imag(aElli),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',2);    hold off;
    hold on;   plot(aTime(1:N-50),real(aElli(1:N-50)),'Color',gwlGetColor(0),'LineWidth',2);    hold off;
    line([min(aTime),max(aTime)],[0,0],'Color','black');
    line([min(aTime),min(aTime)],[-aMax,aMax],'Color','black');
    legend(gwlGetNotation('MSIG','T','x'),gwlGetNotation('MSIG','T','z'));
    gwlText(max(aTime)*0.91,-0.1,gwlGetNotation('TIME'));
    aInd1 = 130; line([aTime(aInd1),aTime(aInd1)],[2,9],'Color','black');
    aInd2 = 247; line([aTime(aInd2),aTime(aInd2)],[2,9],'Color','black');
    line([aTime(aInd1),aTime(aInd2)],[3,3],'Color','black');
    gwlText(0.47*(aTime(aInd1)+aTime(aInd2)),2.5,gwlGetNotation('EPAR','PDIFF'));
    set(gca,'Visible','Off');

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc exmETparFig1;
    
