function createWavelets()
% createWavelets(): Reading and plotting wavelets representation, that is calculated using gwlWavelets porogram

path(path, '../../mshell');

%---------------------------------------------------------------------------
figure(1)
aScale = 2;
[aTime,aWavRe,aWavIm,aFreq,aFourRe,aFourIm]=textread('waveletCauchy.dat','%f %f %f %f %f %f');
aWav = (aWavRe+i*aWavIm)/aScale;
aFour = (aFourRe+i*aFourIm)/aScale;
gwlPlotFunction(aTime,real(aWav),0.07,0.3,0.4,0.3,-3,3,-1,1,gwlGetNotation('TIME'),gwlGetNotation('GSIG','T'),'(a)');
    hold on; plot(aTime,imag(aWav),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend(gwlGetNotation('GSIG','RET'),gwlGetNotation('GSIG','IMT'))
gwlPlotFunction(aFreq,abs(aFour),0.55,0.3,0.4,0.3,-5,5,0,max(abs(aFour)),gwlGetNotation('FREQ'),gwlGetNotation('GSIG','F'),'(b)');

%---------------------------------------------------------------------------
figure(2)
aScale = 0.5;
[aTime,aWavRe,aFreq,aFourRe,aFourIm]=textread('waveletHaar.dat','%f %f %f %f %f');
aWav = (aWavRe+i*aWavIm)/aScale;
aFour = (aFourRe+i*aFourIm)/aScale;
gwlPlotFunction(aTime,real(aWav),0.07,0.3,0.4,0.3,-5,5,-2,2,gwlGetNotation('TIME'),gwlGetNotation('GSIG','T'),'(a)');
gwlPlotFunction(aFreq,imag(aFour),0.55,0.3,0.4,0.3,-5,5,-max(abs(aFour)),max(abs(aFour)),gwlGetNotation('FREQ'),gwlGetNotation('GSIG','F'),'(b)');

%---------------------------------------------------------------------------
figure(3)
aScale = 2;
[aTime,aWavRe,aWavIm,aFreq,aFourRe,aFourIm]=textread('waveletMorlet.dat','%f %f %f %f %f %f');
aWav = (aWavRe+i*aWavIm)/aScale;
aFour = (aFourRe+i*aFourIm)/aScale;
gwlPlotFunction(aTime,real(aWav),0.07,0.3,0.4,0.3,-3,3,-1,1,gwlGetNotation('TIME'),gwlGetNotation('GSIG','T'),'(a)');
    hold on; plot(aTime,imag(aWav),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend(gwlGetNotation('GSIG','RET'),gwlGetNotation('GSIG','IMT'))
gwlPlotFunction(aFreq,abs(aFour),0.55,0.3,0.4,0.3,-5,5,0,max(abs(aFour)),gwlGetNotation('FREQ'),gwlGetNotation('GSIG','F'),'(b)');

%---------------------------------------------------------------------------
figure(4)
gwlPlotFunction(aTime*aScale,real(aWav),0.07,0.3,0.4,0.3,-10,10,-1,1,gwlGetNotation('TIME'),gwlGetNotation('MSIG','T'),'(a)');
    hold on;  plot(aTime*aScale,abs(aWav),'--black','LineWidth',1);  hold off;
    hold on;  plot(aTime*aScale,-abs(aWav),'--black','LineWidth',1);  hold off;
gwlPlotFunction(aFreq/aScale,abs(aFour),0.55,0.3,0.4,0.3,0,5,0,max(abs(aFour)),gwlGetNotation('CFREQ'),gwlGetNotation('MSIG','CF',0),'(b)');

%---------------------------------------------------------------------------
figure(5)
aScale = 2;
[aTime,aWavRe,aFreq,aFourRe,aFourIm]=textread('waveletReMorlet.dat','%f %f %f %f %f');
aWav = (aWavRe+i*aWavIm)/aScale;
aFour = (aFourRe+i*aFourIm)/aScale;
gwlPlotFunction(aTime,real(aWav),0.07,0.3,0.4,0.3,-5,5,-1,1,gwlGetNotation('TIME'),gwlGetNotation('GSIG','T'),'(a)');
gwlPlotFunction(aFreq,abs(aFour),0.55,0.3,0.4,0.3,-5,5,0,max(abs(aFour)),gwlGetNotation('FREQ'),gwlGetNotation('GSIG','F'),'(b)');

%---------------------------------------------------------------------------
figure(6)
[aTime,aWavRe,aWavIm,aFreq,aFourRe,aFourIm]=textread('waveletShanon.dat','%f %f %f %f %f %f');
aWav = (aWavRe+i*aWavIm);
aFour = (aFourRe+i*aFourIm);
gwlPlotFunction(aTime,real(aWav),0.07,0.3,0.4,0.3,-4,4,-1,1,gwlGetNotation('TIME'),gwlGetNotation('GSIG','T'),'(a)');
    hold on; plot(aTime,imag(aWav),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1); hold off;
    legend(gwlGetNotation('GSIG','RET'),gwlGetNotation('GSIG','IMT'))
gwlPlotFunction(aFreq,abs(aFour),0.55,0.3,0.4,0.3,-5,5,0,max(abs(aFour))*1.1,gwlGetNotation('FREQ'),gwlGetNotation('GSIG','F'),'(b)');

%---------------------------------------------------------------------------
pause(0.00001);
clear all;

print -f1 -r600 -depsc waveletCauchy;
print -f2 -r600 -depsc waveletHaar;
print -f3 -r600 -depsc waveletMorletFig1;
print -f4 -r600 -depsc waveletMorletFig2;
print -f5 -r600 -depsc waveletReMorlet;
print -f6 -r600 -depsc waveletShanon;



