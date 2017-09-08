function SynthSigB()
% In this example we compare of the dispersion models formulated in the Fourier 
% and wavelet spaces for complex synthetic signal. As an approximation of the 
% phase velocity function, we use the three-parameter exponential approximation.
% 
% [1] M.A.Kulesh, M.S.Diallo and M.Holschneider Wavelet analysis of ellipticity, 
%     dispersion, and dissipation properties of Rayleigh waves // Acoustical Physics. 
%     V. 51. No. 4. P. 421-434 (2005).

%---------------------------------------------------------------------------
path(path, '../../mshell');
aFreqName = 'freq.dat';
aModelName = 'model.dat';
aSignalName = 'signal.dat';
aSignalPropName = 'signalprop.dat';
aSpectrName = 'spectrum.dat';
aYmax = 0.05;

%---------------------------------------------------------------------------
aFreq = gwlCreateAxis(256,0.001,13,'lin',aFreqName,'Frequency');
[aFreq, aModel] = gwlDispModel(aFreqName, 'vel', '1300,200,6', 'polin', '0,0.00003',aModelName);
[aTime, aSignal, aParSig] = gwlSignalRead(2,'SynthSigB.asc','func','--istime',aSignalName,'Synthetic complex signal');

%---------------------------------------------------------------------------
figure(1);
gwlExec('gwlDiffeoDisp',[' --infile=' aSignalName ' --outfile=' aSignalPropName ' --model=' aModelName ' --step=1 --dist=3500']);
fid = fopen(aSignalPropName,'r'); [aTimeProp1,aSignalProp1]=gwlReadSignal(fid); fclose(fid);
gwlCreateAxis(256,0.001,13,'lin --sign=full',aFreqName,'Frequency');
gwlCwt(2, aSignalName, aFreqName, 1, 'morlet', 2, aSpectrName,'wavelet spectrum before diffeomorphism');
gwlExec('gwlDiffeoDisp',[' --infile=' aSpectrName ' --outfile=' aSpectrName ' --model=' aModelName ' --prop=1 --step=1 --dist=3500']);
[aTimeProp2,aSignalProp2,aParSig] = gwlIwt(2, aSpectrName, 'delta');
[aFreqFT,aFour1] = gwlFourTrans('mat',aSignalProp2(:,1),2,aParSig.aSample);
[aFreqFT,aFour2] = gwlFourTrans('mat',aSignalProp2(:,2),2,aParSig.aSample);

gwlPlotFunction(aFreq,aModel(:,3),0.13,0.5,0.35,0.23,-max(aFreq),max(aFreq),1200,1500,gwlGetNotation('FREQ'),gwlGetNotation('DISP','VEL'),'(a)');
    hold on;    plot(aFreq,aModel(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    hold on;    plot(-aFreq,aModel(:,3),'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);    hold off;
    hold on;    plot(-aFreq,aModel(:,4),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    legend(gwlGetNotation('DISP','CP','F'),gwlGetNotation('DISP','CG','F'));
gwlPlotFunction(aFreqFT,abs(aFour1),0.56,0.5,0.35,0.23,-max(aFreq),max(aFreq),0,max(abs(aFour1)),gwlGetNotation('FREQ'),'','(b)');
    hold on;    plot(aFreqFT,abs(aFour2),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;
gwlPlotFunction(aTimeProp1,real(aSignalProp1(:,1)),0.07,0.28,0.9,0.13,min(aTimeProp1),max(aTimeProp1),-aYmax,aYmax,'',gwlGetNotation('CSIG','T',1),'(c)');
    hold on;    plot(aTimeProp1,imag(aSignalProp1(:,1)),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
gwlPlotFunction(aTimeProp1,real(aSignalProp1(:,2)),0.07,0.15,0.9,0.13,min(aTimeProp1),max(aTimeProp1),-aYmax,aYmax,gwlGetNotation('TIME'),gwlGetNotation('CSIG','T',2),'(d)');
    hold on;    plot(aTimeProp1,imag(aSignalProp1(:,2)),'Color',gwlGetColor(0),'LineStyle','--','LineWidth',1);    hold off;
    hold on;    plot(aTimeProp2,real(aSignalProp2(:,2)),'Color',gwlGetColor(1),'LineStyle','-','LineWidth',1);    hold off;
    hold on;    plot(aTimeProp2,imag(aSignalProp2(:,2)),'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);    hold off;

%---------------------------------------------------------------------------
figure(2);
[aTimeWT,aFreqWT,aCwtAbs] = gwlConvert('3','',aSpectrName);
[aTimeWT,aFreqWT,aCwtArg] = gwlConvert('5','--filter=10',aSpectrName);
[ASizex,ASizey,ASizez]=size(aCwtAbs);
for i=1:ASizex 
    for j=1:ASizey
        if(j<(ASizey/2-5)) 
            AWavarg(i,j)=aCwtArg(i,j,1); 
            AWavabs(i,j)=aCwtAbs(i,j,1); 
        else if(j>(ASizey/2+5)) 
                 AWavarg(i,j)=aCwtArg(i,j,2); 
                 AWavabs(i,j)=aCwtAbs(i,j,2); 
             else
                 AWavarg(i,j)=-pi; 
                 AWavabs(i,j)=0; 
             end;
        end;
    end;
end;

gwlPlotImage(aTimeWT,aFreqWT,AWavabs,0.07,0.42,0.9,0.3,'',gwlGetNotation('FREQ'),'(e)');
    line([0,max(aTimeWT)],[0,0],'Color','black');
    gwlText(1.4,11,gwlGetNotation('CSIG','WABS',1));
    gwlText(4.3,11,gwlGetNotation('CSIG','WABS',2));
gwlPlotImage(aTimeWT,aFreqWT,AWavarg,0.07,0.1,0.9,0.3,gwlGetNotation('TIME'),gwlGetNotation('FREQ'),'(f)');
    line([0,max(aTimeWT)],[0,0],'Color','black');
    gwlText(1.4,11,gwlGetNotation('CSIG','WARG',1));
    gwlText(4.3,11,gwlGetNotation('CSIG','WARG',2));

  
%---------------------------------------------------------------------------
pause(0.00001);
delete(aFreqName);  delete(aModelName);  delete(aSignalName);  delete(aSignalPropName);  delete(aSpectrName); 
clear all;

print -f1 -r600 -depsc SynthSigBFig1;
print -f2 -r600 -depsc SynthSigBFig2;
