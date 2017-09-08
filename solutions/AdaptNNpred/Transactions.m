function Transactions() 
% In this example, the fragments of dirty bank data are presented. This data is 
% related to the processing of bank-to-bank transactions and commercial 
% transactions, which occur during week days. Calculations for both series are 
% done for 10% of series length.
%
% [1] K. Kurennaya, M. Kulesh, M. Holschneider Adaptive metrics in the nearest 
%     neighbours method // Preprint Series DFG SPP 1114, University of Bremen. 
%     Preprint 139 (2006).
% [2] M. Kulesh, et al., Adaptive metrics in the nearest neighbours 
%     method, Physica D (2007), doi:10.1016/j.physd.2007.08.019
%
% FIGURE 1. (a) Real data -- payments based on paper-documents, vertical line 
% shows the prediction start. (b) Zoom of predicted values
%
% FIGURE 2. (a) Real data -- number of transactions, vertical line shows the 
% prediction start. (b) Zoom of predicted values

%----------------------------------------------------------------------------
path(path, '../../mshell');

%----------------------------------------------------------------------------
figure(1); 
gwlExec('gwlNNpred','--infile=TransactionsSig1.asc --outfile=TransactionsSig1.dat --outtype=1 --slength=500 --plength=50 --l3 --dim=20 --nncount=2 --pcoeff=1 --name="Transactions 1"');
localPlotOneSerie('TransactionsSig1.asc', 'TransactionsSig1.dat',2); 

%----------------------------------------------------------------------------
figure(2); 
gwlExec('gwlNNpred','--infile=TransactionsSig2.asc --outfile=TransactionsSig2.dat --outtype=1 --slength=350 --plength=34 --l3 --dim=40 --nncount=2 --pcoeff=1 --name="Transactions 2"');
localPlotOneSerie('TransactionsSig2.asc', 'TransactionsSig2.dat',0.5); 

%---------------------------------------------------------------------------
pause(0.00001);
print -f1 -r600 -depsc TransactionsFig1;
print -f2 -r600 -depsc TransactionsFig2;

%---------------------------------------------------------------------------
% Local functions
%---------------------------------------------------------------------------
function localPlotOneSerie(aSourceFile, aResFile, aMin) 
aSourceFile
aSource = load(aSourceFile,'-ascii');
aRes    = load(aResFile,'-ascii'); 
aStart  = length(aSource)-length(aRes)+1;
aEnd    = length(aSource);
calcMAPE(aSource, aRes, aStart, aEnd)

aMax = max(aSource);
aAmpl=aMin:(aMax-aMin)/100:aMax;
gwlPlotFunction(1:aEnd,aSource,0.07,0.57,0.9,0.4,1,aEnd,aMin,aMax,'Time','Data series','(a)');
    hold on; 
    plot(aStart:aEnd,aRes,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);
    hold off;
    line([aStart,aStart],[aMin,aMax],'Color',gwlGetColor(1));
    grid off;
aMin = min(aSource(aStart:aEnd));
aMax = max(aSource(aStart:aEnd));
gwlPlotFunction(aStart:aEnd,aSource(aStart:aEnd),0.07,0.07,0.9,0.4,aStart,aEnd,aMin,aMax*1.1,'Time','Data series','(b)');
    hold on; 
    plot(aStart:aEnd,aRes,'Color',gwlGetColor(1),'LineStyle','--','LineWidth',1);
    hold off;
    legend('Source data','Prediction using L_A');
    grid off;

