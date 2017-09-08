function RealSigA()
% We apply the propagation operator to a real data consisting of seismic data 
% obtained from a shallow seismic experiment. From the recorded seismograms we 
% selected some traces on successive stations for our test. The Hammer strike is 
% used as seismic source with equidistant geophones (distance between stations 2 
% m) as receiver. First we determine the phase velocity, the group velocity and 
% the attenuation coefficient using Fourier-based methods. Instead of using an 
% approach based on two stations, we considered the parameterization of the 
% wavenumber and attenuation functions and estimated them via a properly defined 
% cost function that involves all the preselected traces. 
% 
% Taking the traces closest to the source for reference, we then try to reproduce 
% the observed signal at further stations, using the propagation operator with the 
% attenuation and the wavenumber (and its derivative) estimated above.
% 
% [1] M.Kulesh, M.Holschneider, M.S.Diallo, Q.Xie and F.Scherbaum Modeling of wave 
%     dispersion using continuous wavelet transforms // Pure and Applied Geophysics. 
%     V. 162. No. 5. P. 843-855 (2005).
% 
% FIGURE 1. Test of the wavelet propagator on experimental data. (a) The 
% wavenumber curve, (b) the phase and group velocity and (c) the attenuation 
% coefficient are determined within the Fourier domain by a minimization process 
% that seeks to find the parameterized wavenumber and attenuation which minimize a 
% properly chosen objective function. (d) Comparison between the original traces 
% (red) with those reproduced by propagating the reference trace (bottom trace on 
% the panel) to the successive receiver stations (black lines). Note the fairly 
% good agreement between the original and reproduced signal especially for the 
% strong arrival associated with the surface waves.
% 
% FIGURE 2. (a) The modulus of the wavelet transforms of reference waveforms 
% (bottom Fig. 1(d) and (b) its phase image. (c) The modulus of the wavelet 
% transforms of last waveforms (fifth waveform in Figures 1(d) and (d) its
% phase image for the frequency range 20 to 60 Hz.
    
%---------------------------------------------------------------------------
path(path, '../../mshell');
aFreqName = 'freq.dat';
aModelName = 'model.dat';
aSignalName = 'signal.dat';
aSpectrName = 'spectrum.dat';

%---------------------------------------------------------------------------
aFreq = gwlCreateAxis(128,20,60,'lin',aFreqName,'Frequency');
[aFreq, aModel] = gwlDispModel(aFreqName, 'polin', '0.0,-1.6392e-02,+3.9406e-03,-1.8138e-04,+2.4596e-06,+5.7903e-09,-1.2064e-10', 'polin', '+2.8130e+00,-4.9107e-01,+6.3838e-02,-2.87607e-03,+5.2688e-05,-3.1584e-07',aModelName);
[aTime, aSignal] = gwlSignalRead(1,'RealSigA.asc','seis','--chan=17 --smplfreq=4000 --to2p',aSignalName,'Experimental signal');
[aTimeProp1, aSignalProp1] = gwlSignalRead(1,'RealSigA.asc','seis','--chan=17,18,19,20,21,22 --smplfreq=4000 --to2p');
gwlCwt(1, aSignalName, aFreqName, 1, 'morlet', 1, aSpectrName,'wavelet spectrum before diffeomorphism');
gwlExec('gwlDiffeoDisp',[' --infile=' aSpectrName ' --outfile=' aSpectrName ' --model=' aModelName ' --step=5 --dist=2']);
[aTimeProp2,aSignalProp2] = gwlIwt(1, aSpectrName, 'delta');
[aTime,aFreq,aCwtAbs] = gwlConvert('3','',aSpectrName);
[aTime,aFreq,aCwtArg] = gwlConvert('5','--filter=20',aSpectrName);

%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aFreq,aModel(:,1),0.07,0.7,0.26,0.25,min(aFreq),max(aFreq),0,1,gwlGetNotation('FREQ'),gwlGetNotation('DISP','WN','F'),'(a)');
gwlPlotFunction(aFreq,aModel(:,3),0.39,0.7,0.26,0.25,min(aFreq),max(aFreq),0,300,gwlGetNotation('FREQ'),gwlGetNotation('DISP','VEL'),'(b)');
    hold on;    plot(aFreq,aModel(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend(gwlGetNotation('DISP','CP','F'),gwlGetNotation('DISP','CG','F'));
gwlPlotFunction(aFreq,aModel(:,5),0.71,0.7,0.26,0.25,min(aFreq),max(aFreq),0,2,gwlGetNotation('FREQ'),gwlGetNotation('DISP','ATN','F'),'(c)');
gwlPlotSeis(aTimeProp1,aSignalProp1,0.07,0.1,0.9,0.5,min(aTimeProp1),max(aTimeProp1),1,gwlGetNotation('TIME'),gwlGetNotation('TRN'),0,'(d)');
    gwlPlotSeisAdd(aTimeProp2,aSignalProp2,1,0,1);
    grid off;

%---------------------------------------------------------------------------
figure(2);
gwlPlotImage(aTime,aFreq,aCwtAbs(:,:,1),0.07,0.75,0.9,0.205,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('SIG','WABS','1')]);
gwlPlotImage(aTime,aFreq,aCwtArg(:,:,1),0.07,0.53,0.9,0.205,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('SIG','WARG','1')]);
gwlPlotImage(aTime,aFreq,aCwtAbs(:,:,5),0.07,0.31,0.9,0.205,'',gwlGetNotation('FREQ'),['(c) ' gwlGetNotation('SIG','WABS','5')]);
gwlPlotImage(aTime,aFreq,aCwtArg(:,:,5),0.07,0.09,0.9,0.205,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),['(d) ' gwlGetNotation('SIG','WARG','5')]);

%---------------------------------------------------------------------------
pause(0.00001);
delete(aFreqName);  delete(aModelName);  delete(aSignalName);  delete(aSpectrName);
clear all;

print -f1 -r600 -depsc RealSigAFig1;
print -f2 -r600 -depsc RealSigAFig2;
