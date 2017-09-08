function exmCwtHER()
% exmCwtHER(): Two component polarization filter of a geomagnetic record.
%
% [1] M. Kulesh, M. Nose and M. Holschneider. Polarization Analysis of Pi2 
%     Pulsations Using Continuous Wavelet Transform // Eos Trans. AGU, 87(52), Fall 
%     Meet. Suppl., Abstract SM43D-02 (2006).
%
% FIGURE 1. We get here two orthogonal components of geomagnetic record from 
% Hermanius station. From these components we construct a complex signal "Z(t)". 
% Next, we transform this signal to time-frequency domain. Finally, we plot 
% separately the image of modulus and the image of complex phases of wavelet 
% coefficients. The source signal is complex; therefore, we have the wavelet 
% spectrum for positive as well as for negative frequencies. We call it 
% progressive and regressive spectrum. 
% 
% FIGURE 2. After the continuous wavelet transform, we have a complex wavelet 
% spectrum WgZ(t,f). We delete from this spectrum the points or coefficients, which 
% do not correspond to the ratio and the tilt angle permissible by the filtering 
% intervals "P_rho" and "P_theta". After this, we have a new wavelet spectrum F_e.

%---------------------------------------------------------------------------
path(path, '../../mshell');
aFreqName = 'freq.dat';
gwlCreateAxis(128,0.0005,0.04,'lin --sign=full',aFreqName,'Frequency');
aCat = 40;
aX = 0.07;
aY = 0.9;

%---------------------------------------------------------------------------
figure(1); 
evalOneStation('1997-01-14', 'HER', '0,1', [42300 44348], aFreqName, 1);
fid = fopen('1997-01-14HERsig.dat','r');    [aTime,aSig] = gwlReadSignal(fid);    fclose(fid);  
[aTime,aFreq,aCwt] = gwlConvert('3,5','--filter=5','1997-01-14HERcwt.dat');
aTime = aTime./3600;
aWgArg = 180.0*aCwt(:,:,2)/pi;
gwlPlotFunction(aTime,real(aSig),aX,aY-0.1,0.77,0.08,min(aTime),max(aTime),min(real(aSig)),max(real(aSig)),'',gwlGetNotation('MSIG','T','n'),'(a)');
gwlPlotFunction(aTime,imag(aSig),aX,aY-0.185,0.77,0.08,min(aTime),max(aTime),min(imag(aSig)),max(imag(aSig)),'',gwlGetNotation('MSIG','T','e'));
gwlPlotImage(aTime, aFreq, aCwt(:,:,1), aX,aY-0.397,0.9,0.2, '', gwlGetNotation('FREQ'), ['(b) ' gwlGetNotation('CSIG','WABS')]);
   clb = colorbar;
   set(clb,'FontSize',7);
gwlPlotImage(aTime, aFreq, aWgArg, aX,aY-0.61,0.9,0.2, gwlGetNotation('TIME','hours'), gwlGetNotation('FREQ'), ['(c) ' gwlGetNotation('CSIG','WARG')]);
   clb = colorbar;
   set(clb,'FontSize',7);

%---------------------------------------------------------------------------
figure(2); 
gwlExec('gwlET2DFilter',[' --infile=1997-01-14HERcwt.dat --outfile=1997-01-14HERcwtfil.dat --filter=linhor,0.4,0.25 --elli=1997-01-14HERpar.dat --name="filtered cpectrum"']);
[aTime,aFreq,aCwtFil] = gwlConvert('3,5','--filter=5','1997-01-14HERcwtfil.dat');
aTime = aTime./3600;
gwlPlotFunction(aTime,real(aSig),aX,aY-0.1,0.77,0.08,min(aTime),max(aTime),min(real(aSig)),max(real(aSig)),'',gwlGetNotation('MSIG','T','n'),'(a)');
gwlPlotFunction(aTime,imag(aSig),aX,aY-0.185,0.77,0.08,min(aTime),max(aTime),min(imag(aSig)),max(imag(aSig)),'',gwlGetNotation('MSIG','T','e'));
gwlPlotImage(aTime, aFreq, aCwt(:,:,1), aX,aY-0.397,0.9,0.2, '', gwlGetNotation('FREQ'), ['(b) ' gwlGetNotation('CSIG','WABS')]);
   clb = colorbar;
   set(clb,'FontSize',7);
gwlPlotImage(aTime, aFreq, aCwtFil(:,:,1), aX,aY-0.61,0.9,0.2, gwlGetNotation('TIME','hours'), gwlGetNotation('FREQ'), ['(c) F_e' gwlGetNotation('CSIG','WABS')]);
   clb = colorbar;
   set(clb,'FontSize',7);

%---------------------------------------------------------------------------
pause(0.00001);
delete('*.dat');
clear all;

print -f1 -r600 -depsc exmCwtHERFig1;
print -f2 -r600 -depsc exmCwtHERFig2;
