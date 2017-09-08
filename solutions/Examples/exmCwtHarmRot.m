function exmCwtHarmRot()
% These pictures was used during the talk at conference:
% Diallo M.S., Kurennaya K., Kulesh M. and Holschneider M
% Elliptic properties of surface elastic waves in wavelet domain 
% (in Russian) // Book of abstracts of XIV Winter School on 
% Continuous Media Mechanics (28 February - 3 March 2005, Perm). P. 100.

%---------------------------------------------------------------------------
path(path, '../../mshell');

aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';
aYmax = 10.0;
 
gwlCreateAxis(1024,0,5.115,'lin',aTimeName,'Time');
gwlCreateAxis(256,0.0001,15,'lin --sign=full',aFreqName,'Frequency');

[aTime,aSignal,aParSig] = gwlSignalGen(2,aTimeName,'harmrot','2.0,7.0,1.0,2.0,5.0,5.0,0.318',aSignalName,'Rotated harmonic function');
gwlPlotFunction(aTime,real(aSignal),0.07,0.8,0.9,0.15,min(aTime),max(aTime),-aYmax,aYmax,'',gwlGetNotation('CSIG','T'),'(a)');
    hold on; plot(aTime,imag(aSignal),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1); hold off;

gwlCwt(2, aSignalName, aFreqName, 2, 'cauchy', 10, aSpectrName);
[aTime,aFreq,aCwt] = gwlConvert('3,5','--filter=10',aSpectrName);

gwlPlotImage(aTime, aFreq, aCwt(:,:,1), 0.07,0.64,0.9,0.15, '', gwlGetNotation('FREQ'), '(b)');
      line([0,max(aTime)],[0,0],'Color','black');
gwlPlotImage(aTime, aFreq, aCwt(:,:,2), 0.07,0.48,0.9,0.15, gwlGetNotation('TIME'), gwlGetNotation('FREQ'), '(c)');
      line([0,max(aTime)],[0,0],'Color','black');

gwlPlotFunction(real(aSignal),imag(aSignal),0.07,0.20,0.42,0.2,-aYmax,aYmax,-aYmax,aYmax,gwlGetNotation('CSIG','RET'),gwlGetNotation('CSIG','IMT'),'(d)');

[aFreqFT,aFour] = gwlFourTrans('mat',aSignal,2,aParSig.aSample);
gwlPlotFunction(aFreqFT,abs(aFour),0.55,0.20,0.42,0.2,min(aFreq),max(aFreq),0,max(abs(aFour)),gwlGetNotation('FREQ'),gwlGetNotation('CSIG','FM'),'(e)');
      
%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete(aSignalName); delete(aSpectrName); 
clear all;

print -f1 -r600 -depsc exmCwtHarmRot;
