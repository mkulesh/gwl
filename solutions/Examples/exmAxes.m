function exmAxes()
% exmAxes(): An example of linear and logarithimic axes

path(path, '../../mshell');

[aAxLin,aParLin] = gwlCreateAxis(1000,1,100,'lin --sign=full','axislin.dat','Linear axis');
aParLin
gwlPlotFunction(1:length(aAxLin), aAxLin, 0.07,0.3,0.4,0.4, 1,length(aAxLin),min(aAxLin),max(aAxLin),gwlGetNotation('IND'),gwlGetNotation('AVAL'),'(a)',aParLin.aName);

[aAxLog,aParLog] = gwlCreateAxis(1000,1,100,'log --sign=full','axislog.dat','Logarithmic axis');
aParLog
gwlPlotFunction(1:length(aAxLog), aAxLog, 0.55,0.3,0.4,0.4, 1,length(aAxLog),min(aAxLog),max(aAxLog),gwlGetNotation('IND'),gwlGetNotation('AVAL'),'(b)',aParLog.aName);

%---------------------------------------------------------------------------
pause(0.00001);
delete('axislin.dat'); delete('axislog.dat'); 
clear all;
