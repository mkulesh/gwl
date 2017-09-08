function aAns = calcMAPE(aSource, aRes, aStart, aEnd) 

aDeltaSumm=0;
for k=aStart:aEnd
    aDeltaSumm = aDeltaSumm+abs((aSource(k)-aRes(k-aStart+1))/aSource(k));
end;    
aAns.length = 100*length(aRes)/(length(aSource)-length(aRes));
aAns.MAPE = 100*aDeltaSumm/(aEnd-aStart);
