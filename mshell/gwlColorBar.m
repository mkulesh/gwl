function aColorBar = gwlColorBar
% aColorBar = gwlColorBar
% This procedure returns a link to color bar object and can be used after gwlPlotImage method. 
% The way to get such a link is different in a Windows and in a Linux and Matlab versions.
%
% Examples: 
% To modify the YTick's of color bar, the following code can be used:
% gwlPlotImage(...);
%    colorbar;
%    set(gwlColorBar,'YTick',0:0.2:1);
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

str = computer;
if strcmp(str,'PCWIN') == 1
    obj=get(gcf,'Children');
    aColorBar = obj(2);
else
    obj=get(gcf,'Children');
    aColorBar = obj(1);
end;
    