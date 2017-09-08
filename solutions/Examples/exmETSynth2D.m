function exmETSynth2D()

%---------------------------------------------------------------------------
path(path, '../../mshell');

aTimeName = 'time.dat';
aFreqName = 'freq.dat';
aSignalName = 'signal.dat';
aSpectrName = 'cwt.dat';
aElliparName = 'elli.dat';
aFilteredName = 'filtered.dat';
aYmax = 10.0;
aTmin = 0.0;
aTmax = 5.0;
aImax = 1001;
aAlg = 'comple'; 

%---------------------------------------------------------------------------
gwlCreateAxis(1024,0,5.115,'lin',aTimeName,'Time');
[aTime,aSignal] = prepareSignal(aTimeName, aSignalName);

if(aAlg == 'acovar') 
    gwlCreateAxis(128,0.01,30,'lin',aFreqName,'Frequency');
    gwlExec('gwlConvert',[' --infile=' aSignalName ' --outfile=' aSignalName ' --comp=1,2']);
    gwlCwt(1, aSignalName, aFreqName, 0, 'cauchy', 5, aSpectrName);
    gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=acovar --filter=10 --name="polarization properties"']);
    gwlExec('gwlET2DFilter',[' --infile=' aSpectrName ' --outfile=' aFilteredName ' --type=acovar --filter=ellivert,0.785,0.1 --name="filtered spectrum"']);
    [aTimeInv,aSignalInv] = gwlIwt(1, aFilteredName, 'delta');
    aSignalInv = aSignalInv(:,1)+i*aSignalInv(:,2);
else
    gwlCreateAxis(128,0.01,30,'lin --sign=full',aFreqName,'Frequency');
    gwlCwt(2, aSignalName, aFreqName, 0, 'cauchy', 5, aSpectrName);
    gwlExec('gwlET2D',[' --infile=' aSpectrName ' --outfile=' aElliparName ' --type=complex --filter=15 --name="polarization properties"']);
    gwlExec('gwlET2DFilter',[' --infile=' aSpectrName ' --outfile=' aFilteredName ' --type=complex --filter=ellivert,0.785,0.1 --name="filtered spectrum"']);
    [aTimeInv,aSignalInv] = gwlIwt(2, aFilteredName, 'delta');
    % Alternative filtering method
    % gwlExec('gwlET2DFilter',[' --infile=' aSignalName ' --outfile=' aFilteredName ' --type=complex --filter=ellivert,0.785,0.1 --name="filtered spectrum" --wavelet=cauchy --wavpar=5 --freq=' aFreqName]);
    % fid = fopen(aFilteredName,'r'); 
    % [aTimeInv,aSignalInv]=gwlReadSignal(fid); 
    % fclose(fid);
end;

[aTime,aFreq,aCwt] = gwlConvert('3','',aSpectrName);
[aTimeE,aFreqE,aCwtE] = gwlConvert('1,2,4,14','',aElliparName);
[aTime2,aRmax2,aRmin2,aPhidiff2,aTilt2,aRatio2] = calcElliSignal(aSignalName,2);
[aTime3,aRmax3,aRmin3,aPhidiff3,aTilt3,aRatio3] = calcElliSignal(aElliparName,3);
 
%---------------------------------------------------------------------------
figure(1);
gwlPlotFunction(aTime,real(aSignal),0.07,0.82,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','x'),'(a)');
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotFunction(aTime,imag(aSignal),0.07,0.68,0.9,0.14,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','z'),'');
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotImage(aTime(1:aImax), aFreq, aCwt(:,1:aImax,1), 0.07,0.43,0.9,0.24, gwlGetNotation('TIME'), gwlGetNotation('FREQ'), '(b)');
    line([0,aTmax],[0,0],'Color','black');
    if(aAlg=='acovar') 
        gwlText(0.5,15,gwlGetNotation('MSIG','WABS','x'));
    else
        gwlText(0.5,15,gwlGetNotation('CSIG','WABSP'));
        gwlText(0.5,-15,gwlGetNotation('CSIG','WABSM'));
    end;
gwlPlotFunction(real(aSignal),imag(aSignal),0.4,0.07,0.24,0.29,-aYmax,aYmax,-aYmax,aYmax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('MSIG','T','z'),'(c)');

%---------------------------------------------------------------------------
figure(2);
aSemiAxes = cat(1,aCwtE(:,:,1).^2,aCwtE(:,:,2).^2);
gwlPlotImage(aTime(1:aImax), aFreq, aSemiAxes(:,1:aImax), 0.07,0.76,0.9,0.22 ,'',gwlGetNotation('FREQ'),'(a)');
    AF1 = min(aFreq);
    AF2 = max(aFreq);
    set(gca,'YTick',AF1:AF2/3:AF2);
    set(gca,'YTickLabel',{0,AF2/3,2*AF2/3,0,AF2/3,2*AF2/3,AF2});
    gwlText(1.1,20,gwlGetNotation('EPAR','RMIN','WG'));
    gwlText(1.1,-4,gwlGetNotation('EPAR','RMAX','WG'));
    line([0,aTmax],[0,0],'Color','black');
gwlPlotFunction(aTime3,aRatio3,0.07,0.528,0.9,0.22,aTmin,aTmax,0,1,'',gwlGetNotation('EPAR','RATIO','T'),'(b)');
    hold on;   plot(aTime2,aRatio2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    legend('This study','Rene et al. (1986)');
gwlPlotFunction(aTime3,aPhidiff3,0.07,0.296,0.9,0.22,aTmin,aTmax,-pi,pi,'',gwlGetNotation('EPAR','PDIFF','T'),'(c)');
    hold on;   plot(aTime2,aPhidiff2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
    legend('This study','Rene et al. (1986)');
gwlPlotFunction(aTime3,aTilt3,0.07,0.06,0.9,0.22,aTmin,aTmax,-pi,pi,gwlGetNotation('TIME'),gwlGetNotation('EPAR','TILT','T'),'(d)');
    hold on;   plot(aTime2,aTilt2,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-pi:pi/2:pi);
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
    legend('This study','Rene et al. (1986)');

%---------------------------------------------------------------------------
figure(3);
gwlPlotImage(aTimeE(1:aImax), aFreqE, aCwtE(:,1:aImax,3),0.2,0.82,0.72,0.15,'',gwlGetNotation('FREQ'),['(a) ' gwlGetNotation('EPAR','RATIO','WG')]);
    colorbar;
    set(gwlColorBar,'YTick',0:0.2:1);
    set(gwlColorBar,'YTickLabel',{0,0.2,0.4,0.7,0.8,1});
gwlPlotImage(aTimeE(1:aImax), aFreqE, aCwtE(:,1:aImax,4),0.2,0.66,0.72,0.15,'',gwlGetNotation('FREQ'),['(b) ' gwlGetNotation('EPAR','TILT','WG')]);
    colorbar;
    set(gwlColorBar,'YTick',-1.5:0.75:1.5);
    set(gwlColorBar,'YTickLabel',{'-pi/2','-pi/4','0','pi/4','pi/2'});
gwlPlotFunction(aTime,real(aSignal),0.2,0.5,0.614,0.15,aTmin,aTmax,-aYmax,aYmax,'',gwlGetNotation('MSIG','T','x'),'(c)');
    hold on;   plot(aTimeInv,real(aSignalInv),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotFunction(aTime,imag(aSignal),0.2,0.34,0.614,0.15,aTmin,aTmax,-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('MSIG','T','z'),'(d)');
    hold on;   plot(aTimeInv,imag(aSignalInv),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
    set(gca,'YTick',-aYmax:aYmax/2:aYmax);
    set(gca,'YTickLabel',{'',-aYmax/2,0,aYmax/2,''});
gwlPlotFunction(real(aSignal),imag(aSignal),0.41,0.07,0.18,0.21,-aYmax,aYmax,-aYmax,aYmax,gwlGetNotation('MSIG','T','x'),gwlGetNotation('MSIG','T','z'),'(e)');    
    hold on;   plot(real(aSignalInv),imag(aSignalInv),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
colormap hsv;

%---------------------------------------------------------------------------
pause(0.00001);
delete(aTimeName); delete(aFreqName); delete(aSignalName); delete(aSpectrName); delete(aElliparName); delete(aFilteredName); 
clear all;

print -f1 -r600 -depsc exmETSynth2DFig1;
print -f2 -r600 -depsc exmETSynth2DFig2;
print -f3 -r600 -depsc exmETSynth2DFig3;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function [aTime,aRmax,aRmin,aPhidiff,aTilt,aRatio] = calcElliSignal(aSourceName,aType)
aTmpName1 = tempname;
if(aType == 1)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=complex']);
end;
if(aType == 2)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --type=rene']);
end;
if(aType == 3)
    gwlExec('gwlET2D',[' --infile=' aSourceName ' --outfile=' aTmpName1 ' --freq=1,7']);
end;
gwlExec('gwlConvert',[' --infile=' aTmpName1 ' --outfile=' aTmpName1 ' --outtype=1 --comp=1,2,13,14,4 --nomess']);
[aTime,aRmax,aRmin,aPhidiff,aTilt,aRatio] = textread(aTmpName1,'%f %f %f %f %f %f');
delete(aTmpName1);

function [aTime,aSignal] = prepareSignal(aTimeName,aSourceName)
at0 =  5.115/3;
aMode1 = tempname;
aMode2 = tempname;
aMode3 = tempname;
gwlSignalGen(2,aTimeName,'harmon',['1,2,5,1,3,5,1.570796 --mod=0,' num2str(at0) ',1'],aMode1,'Harmonic function');
gwlSignalGen(2,aTimeName,'harmrot',['2,7,1,2,5,5,0.318 --mod=' num2str(at0) ',5.115,1'],aMode2,'Rotated harmonic function');
gwlExec('gwlWavelets',[' --infile=' aTimeName ' --iscmpl --wavelet=morletre --wavpar=1 --time=3.5 --freq=15 --outtype=1 --outfile=' aMode3])
[aTime,aSignal,aPar] = gwlSignalRead(2,aMode3,'func',['--format=ASCII --istime --mult=0.3 --nomess'],aMode3,'Morlet wavelet');
gwlExec('gwlSignalSum',['--infile=' aMode1 ',' aMode2 ',' aMode3 ' --outfile=' aSourceName ' --name="Synthetic signal"']);
fid = fopen(aSourceName,'r');
[aTime,aSignal]=gwlReadSignal(fid);
fclose(fid);
delete(aMode1); delete(aMode2); delete(aMode3);
