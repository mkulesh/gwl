function WorthBasinReal()
% In this example we applied the method to the raw 2-C data from Fort Worth Basin U. S. A.
% Because the data were not yet edited and preprocessed, we cannot carry out detailed 
% interpretation in terms of wave mode as we did for the synthetics. As such, we limit 
% ourselves to a qualitative appraisal of our method in separating the Rayleigh waves 
% from the overall linearly polarized signals - i.e., combined LH & LV. The results for 
% the linearly polarized events are shown in (c) and (d), where wave packets corresponding 
% to the elliptically polarized events (Rayleigh waves among others) have been significantly 
% attenuated. In contrast, Figures (e) and (f) show where wave packets are enhanced. Note the
% improved resolution for the high-frequency event on the vertical component to the right of 
% the receiver at station 80 and between 1 s and 1.5 s in time.
% 
% [1] M.S.Diallo, M.Kulesh, M.Holschneider, F.Scherbaum, F.Adler. 
%     Characterization of polarization  attributes of seismic waves using continuous wavelet 
%     transforms // Geophysics. V. 71. No. 3. P. V67-V77 (2006).
% [2] Diallo M.S., Kulesh M., Holschneider M. and Scherbaum F. Instantaneous polarization 
%     attributes in the time-frequency domain: application to wave field separation // Eos 
%     Trans. AGU, 85(47), Fall Meet. Suppl., Abstract S31B-1063 (2004).
% [3] Diallo M.S., Kurennaya K., Kulesh M. and Holschneider M
%     Elliptic properties of surface elastic waves in wavelet domain (in Russian) // Book 
%     of abstracts of XIV Winter School on Continuous Media Mechanics (28 February - 3 March 
%     2005, Perm). P. 100.
% [4] Mamadou S. Diallo, Michail Kulesh, Matthias Holschneider, Kristina Kurrenaya and Frank 
%     Scherbaum. Estimating polarization attributes with an adaptive covariance method in the 
%     wavelet domain // SEG Technical Program Expanded Abstracts. V. 24. P. 1014-1017 (2005).
%     
% FIGURE 1. Real 2-C shot gather from ForthWorth Basin. (a) The horizontal component. (b) The
% vertical component. (c) and (d) Filtered real shot gathers showing the overall linearly 
% polarized events (combined LH and LV filter). (e) and (f) Filtered real shot gathers showing 
% the overall elliptically polarized events (combined EH and EV). Ratio = 0.15 and tilt = pi/2 
% because no distinction is made between horizontal and vertical polarization.    

%---------------------------------------------------------------------------
path(path, '../../mshell');
aSeisMax = 0.000005;

%---------------------------------------------------------------------------
figure(1);
fid = fopen('BinData/WorthBasinRealsig.dat','r');     [aTime,aSignal]=gwlReadSignal(fid);  fclose(fid);
fid = fopen('BinData/WorthBasinRealfilt(1).dat','r'); [aTime,aLinHor]=gwlReadSignal(fid); fclose(fid);
fid = fopen('BinData/WorthBasinRealfilt(2).dat','r'); [aTime,aEllHor]=gwlReadSignal(fid); fclose(fid);

gwlPlotSeisV(aTime,real(aSignal),0.15,0.69,0.3,0.28,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(a)');
gwlPlotSeisV(aTime,imag(aSignal),0.48,0.69,0.3,0.28,min(aTime),max(aTime),aSeisMax,' ',' ',0,'(b)');

gwlPlotSeisV(aTime,real(aLinHor),0.15,0.38,0.3,0.28,min(aTime),max(aTime),aSeisMax,' ',gwlGetNotation('TIME'),0,'(c)');
gwlPlotSeisV(aTime,imag(aLinHor),0.48,0.38,0.3,0.28,min(aTime),max(aTime),aSeisMax,' ',' ',0,'(d)');
  
gwlPlotSeisV(aTime,real(aEllHor),0.15,0.07,0.3,0.28,min(aTime),max(aTime),aSeisMax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('TIME'),0,'(e)');
gwlPlotSeisV(aTime,imag(aEllHor),0.48,0.07,0.3,0.28,min(aTime),max(aTime),aSeisMax,gwlGetNotation('MSIG','T','z'),' ',0,'(f)');
    

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc WorthBasinReal;
