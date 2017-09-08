function SiahData()

%----------------------------------------------------------------------------
path(path, '../../mshell');

aSourceName = 'SiahData.asc';
aSignalName = 'SiahDataSig.dat';
aSignalInvName = 'SiahDataSigInv.dat';
aFreqName = 'SiahDataFreq.dat';
aVelName = 'SiahDataVel.dat';
aSpectrName = 'SiahDataCwt.dat';
aFKName = 'SiahDataFV.dat';

%----------------------------------------------------------------------------
[aTime,aSeis] = gwlSignalRead(1,aSourceName,'seis',['--format=ASCII --to2p --chan=15,14,13,12,11,10 --smplfreq=5000'],aSignalName,'Section A');
% CWT
gwlCreateAxis(250,80,120,'lin',aFreqName,'Frequency');
gwlCwt(1, aSignalName, aFreqName, 2, 'morlet', 5, aSpectrName);
% CWT filter
[aTimeInv,aSignalInv] = gwlIwt(1, aSpectrName, 'delta', 1, 2, aSignalInvName);
% Signal optimization
gwlCreateAxis(400,0.001,150,'lin','SiahDataFreq1.dat','Frequency');
gwlExec('gwlDispModel', '--infile=SiahDataFreq1.dat --outfile=SiahDataInit.dat --analyt --wn=bspline --wnpar=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 --atn=bspline --atnpar=0.1,0.1,0.1,0.1,0.1 --name="Initial dispersion model"');
gwlExec('gwlOptiSI', '--infile=SiahDataInit.dat --outfile=SiahDataOpt1Par.dat --sig=SiahDataSigInv.dat --osig=SiahDataOpt1Sig.dat --dist=2 --eps=1E-8 --prog --name="Optimization1"');
% Frequency-velocity analysis
gwlCreateAxis(250,50,250,'lin',aVelName,'Velocity');
gwlExec('gwlTransFK', '--infile=SiahDataCwt.dat --outfile=SiahDataFV.dat --vel=SiahDataVel.dat --inter=spline --corr=arg --filter=1 --norm --dist=2 --prog --name="FVimage"');

%----------------------------------------------------------------------------
figure(1);
aAmplMax = 0.05;
gwlPlotSeis(aTime,aSeis,0.07,0.07,0.88,0.9,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),' ',0,'');
    gwlPlotSeisAdd(aTimeInv,aSignalInv,aAmplMax,0,1);
    gwlText(-0.03,40,gwlGetNotation('TRN'),90);
    grid off;

%----------------------------------------------------------------------------
figure(2);
[aTime,aFreq,aCwt] = gwlConvert('5','--filter=1',aSpectrName);
dX = 0.3;
for k=1:3
    gwlPlotImage(aTime, aFreq, aCwt(:,:,k), 0.05+(k-1)*dX,0.55,dX-0.01,0.4, ' ', ' ', '');
    gwlPlotImage(aTime, aFreq, aCwt(:,:,3+k), 0.05+(k-1)*dX,0.1,dX-0.01,0.4, ' ', ' ', '');
end;

%----------------------------------------------------------------------------
figure(3);
aAmplMax = 0.02;
fid = fopen('SiahDataOpt1Sig.dat','r'); [aOpt1Time,aOpt1Sig]=gwlReadSignal(fid); fclose(fid);
gwlPlotSeis(aTimeInv,aSignalInv,0.07,0.07,0.88,0.9,min(aTime),max(aTime),aAmplMax,gwlGetNotation('TIME'),' ',0,'');
    gwlPlotSeisAdd(aOpt1Time,aOpt1Sig,aAmplMax,0,1);
    gwlText(-0.03,40,gwlGetNotation('TRN'),90);
    grid off;

%----------------------------------------------------------------------------
figure(4);
fid = fopen(aFKName,'r'); [aVel,aFreq,aFK]=gwlReadSpectrum(fid); fclose(fid);
    aSlown = fliplr((1./aVel)');
    aFKs = flipud(aFK');
fid = fopen('SiahDataOpt1Par.dat','r'); [aOpt1Freq,aOpt1Par]=gwlReadDispModel(fid); fclose(fid);
    aPar = aOpt1Par(:,3)/(2*pi);
gwlPlotSurface(aFreq,aVel,aFK',2,0,50,0.07,0.07,0.9,0.9,gwlGetNotation('FREQ'),gwlGetNotation('DISP','CP','F'),'');
    line(aOpt1Freq,aPar);

%----------------------------------------------------------------------------
figure(5);
gwlPlotFunction(aOpt1Freq,aPar,0.1,0.1,0.85,0.8,min(aOpt1Freq),max(aOpt1Freq),min(aPar),500,gwlGetNotation('FREQ'),gwlGetNotation('DISP','WN','F'));

%----------------------------------------------------------------------------
print -f1 -r600 -depsc SiahDataFig1;
print -f2 -r600 -depsc SiahDataFig2;
print -f3 -r600 -depsc SiahDataFig3;
print -f4 -r600 -depsc SiahDataFig4;
print -f5 -r600 -depsc SiahDataFig5;
