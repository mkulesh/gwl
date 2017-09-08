function gwlPlotSeisAdd(aTime,aSource, aAmax,aNorm, aColor)
% gwlPlotSeisAdd(aTime,aSource, aAmax,aNorm, aColor)
% This procedure add a multi-channel function to a existing horizontal seismogram
%
% Input parameters:
%   aTime is an array corresponds to the function argument
%   aSource is an array corresponds to the function values
%   aAmax is the maximal amplitude of the function
%   aNorm defines the normalization mode, if aNorm == 1 then all channel will be normalize to 1 
%   aColor defines the color index for new lines (see gwlGetColor for details)
% 
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

[aPointCount,aSigCount]=size(aSource);
for k=1:aSigCount
    aAmpl = aSource(:,k);  
    if aNorm == 1
        maxd = max(aAmpl);
    else
        maxd = 1;
    end;
    hold on;    
    plot(aTime,aAmpl/maxd+aAmax*k,'Color',gwlGetColor(aColor),'LineStyle','--','LineWidth',1);     
    hold off;
end;