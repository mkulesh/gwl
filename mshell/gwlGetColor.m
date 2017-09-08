function aColor = gwlGetColor(aNum)
% aColor = gwlGetColor(aNum)
% This procedure returns a color depends on the its number aNum.
% The color could be from gray-scaled or color palette that is defined
% in gwlGetaColorShema method.
%
% This is an internal procedure of GWL mshell library and is used in plotting
% methods like gwlPlotFunction, gwlPlotImage etc. to make colors of all 
% GWL plots uniform.
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

% gray-scaled palette
if(gwlGetColorShema == 0)
    if(aNum == 0) aColor = [0.0 0.0 0.0];  end;
    if(aNum == 1) aColor = [0.4 0.4 0.4];  end;
    if(aNum == 2) aColor = [0.6 0.6 0.6];  end;
    if(aNum == 3) aColor = [0.5 0.5 0.5];  end;
    if(aNum == 4) aColor = [0.7 0.7 0.7];  end;
end;

% color palette
if(gwlGetColorShema == 1)
    if(aNum == 0) aColor = [1 0 0];        end;
    if(aNum == 1) aColor = [0 0 0];        end;
    if(aNum == 2) aColor = [0 0 1];        end;
    if(aNum == 3) aColor = [0.5 0.5 0.5];  end;
    if(aNum == 4) aColor = [0 1 0];        end;
end;

