function KerpenData()
% KerpenData(): The experimental data set shows in Figure 1 consists of a 2-D shallow seismic 
% survey (stations along a line) at Kerpen, a particular site in the Lower Rhine embayment where
% the buried scarp of a historically active fault is presumed. Presented vertical displacement
% for 48 channels with 2 m inter-receiver spacing were collected using hammer blows as
% seismic source. We selected a seismogram profile with prominent low frequency, high
% amplitude arrivals that correspond to the surface wave arrivals we intend to characterize.
% We selected two subsections from the seismograms for our analysis. These subsections are 
% labeled ‘‘subsection A’’ and ‘‘subsection B’’ in Figure 1. For each subsection, we perform the 
% ‘‘frequency-velocity’’ analysis with the use of correlation between real-valued wavelet phases 
% with a threshold operation  e = 0.1 percent of maximum modulus of wavelet transform.
% See [1] for more details. 
% 
% [1] M. Kulesh, M. Holschneider, M. Ohrnberger, E. Lueck. Modeling of wave dispersion using 
%     continuous wavelet transforms II: wavelet based frequency-velocity analysis // Pure and 
%     Applied Geophysics. Vol. 165. (2008, in press). DOI 10.1007/s00024-008-0299-7
% 
% FIGURE 1. Observed seismograms obtained from a shallow seismic experiment using a sledgehammer 
% as a source. The distance between consecutive stations is 2 m.
% 
% FIGURE 2. The ‘‘frequency-velocity’’ analysis of subsection A.
% 
% FIGURE 3. The ‘‘frequency-velocity’’ analysis of subsection B.
% 
% IMPORTANT! Before run this example, go to ./BinData and run 
% KerpenData.bat from there. The calculation take about 3 hours.

%----------------------------------------------------------------------------
path(path, '../../mshell');

%----------------------------------------------------------------------------
figure(1);
[aTime,aSeis] = gwlSignalRead(1,'BinData/KerpenData.asc','seis',['--format=ASCII --smplfreq=4000']);
gwlPlotSeis(aTime,aSeis,0.07,0.07,0.88,0.9,min(aTime),max(aTime),2.0,gwlGetNotation('TIME'),' ',1,'');
    line([0.49,0.49],[31,41],'Color','black','LineWidth',2);
    gwlText(max(aTime)*1.01,31,'Subsec. A',90);
    line([0.49,0.49],[45,55],'Color','black','LineWidth',2);
    gwlText(max(aTime)*1.01,45,'Subsec. B',90);
    gwlText(-0.03,40,gwlGetNotation('TRN'),90);
    grid off;

%----------------------------------------------------------------------------
figure(2);
Par.SeisName = 'BinData/KerpenDataAsig.dat';
Par.FKName = 'BinData/KerpenDataAarg.dat';
Par.NameMusic = 'MUSIC/music.16-19.xyz';
Par.NameCapon = 'MUSIC/capon.16-19.xyz';
Par.FigLegend = 'Subsec. A';
plotsFKData(Par);

%----------------------------------------------------------------------------
figure(3);
Par.SeisName = 'BinData/KerpenDataBsig.dat';
Par.FKName = 'BinData/KerpenDataBarg.dat';
Par.NameMusic = 'MUSIC/music.23-27.xyz';
Par.NameCapon = 'MUSIC/capon.23-27.xyz';
Par.FigLegend = 'Subsec. B';
plotsFKData(Par);

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc KerpenDataFig1;
print -f2 -r600 -depsc KerpenDataFig2;
print -f3 -r600 -depsc KerpenDataFig3;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function plotsFKData(aPar)
fid = fopen(aPar.SeisName,'r'); [aTime,aSeis]=gwlReadSignal(fid); fclose(fid);
fid = fopen(aPar.FKName,'r'); [aVel,aFreq,aFK]=gwlReadSpectrum(fid); fclose(fid);
[AMusic_x,AMusic_y,AMusic_z]=textread(aPar.NameMusic,'%f %f %f');
[ACapon_x,ACapon_y,ACapon_z]=textread(aPar.NameCapon,'%f %f %f');
FullSize = length(AMusic_x);
YSize = 420;
XSize = FullSize/YSize;
for j=1:XSize
    for i=1:YSize 
        k = (j-1)*YSize + i;
        AMXf(j) = AMusic_x(k);
        AMYf(YSize-i+1) = AMusic_y(k);
        AMZf(YSize-i+1,j) = AMusic_z(k);
        ACZf(YSize-i+1,j) = ACapon_z(k);
    end;
end;
indi1 = 150;
indi2 = 380;
indj1 = 40;
indj2 = 190;
CScale = 2;
for j=indj1:indj2
    for i=indi1:indi2
        AMX(j-indj1+1) = AMXf(j);
        AMY(i-indi1+1) = AMYf(i);
        AMZ(i-indi1+1,j-indj1+1) = AMZf(i,j);
        ACZ(i-indi1+1,j-indj1+1) = ACZf(i,j)/CScale;
    end;
end;
CCount = 8;
aSlown2 = fliplr(1./AMY);
gwlPlotSeis(aTime,aSeis,0.07,0.78,0.9,0.2,min(aTime),max(aTime),2.0,gwlGetNotation('TIME'),gwlGetNotation('TRN'),1,'(a)');
    legend(aPar.FigLegend);
gwlPlotSurface(AMX,aSlown2,flipud((aFK').^20),2,0,50,0.07,0.07,0.44,0.65,gwlGetNotation('FREQ'),['1/' gwlGetNotation('DISP','CP','F')],'(b)');
    hold on;    contour(AMX,aSlown2,flipud(ACZ),CCount);  hold off;
gwlPlotSurface(AMX,aSlown2,flipud(AMZ),2,0,50,0.52,0.07,0.45,0.65,gwlGetNotation('FREQ'),'','(c)');
    hold on;    contour(AMX,aSlown2,flipud(ACZ),CCount);  hold off;

