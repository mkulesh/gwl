function gwlPlotSeisV(aTime,aSource, aX1,aY1,aWidth,aHeight, aTmin,aTmax,aAmax, aXLab,aYLab, aNorm, aText,aTitle)
% gwlPlotSeisV(aTime,aSource, aX1,aY1,aWidth,aHeight, aTmin,aTmax,aAmax, aXLab,aYLab, aNorm [, aText,aTitle])
% This procedure summarize some standard procedures to plot a multi-channel function as a vertical seismogram
%
% Input parameters define the function and its changing interval, the window position of the plot as well as the labels:
%   aTime is an array corresponds to the function argument
%   aSource is an array corresponds to the function values
%   aX1,aY1,aWidth,aHeight describe the window position of the plot (see subplot for details)
%   aTmin,aTmax define the interval of function argument (see set(gca,'Xlim') for details)
%   aAmax is the maximal amplitude of the function
%   aXLab, aYLab are a string variables correspond to the X- and Y- labels of the plot
%   aNorm defines the normalization mode, if aNorm == 1 then all channel will be normalize to 1 
%   aText is a string variable defines the text label in the left-top corner of the plot
%   aTitle is a string variable defines the plot title (see gwlTitle for details)
% 
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aFont = gwlGetFont;
[aPointCount,aSigCount]=size(aSource);
subplot('Position',[aX1,aY1,aWidth,aHeight]);
plot(aTime,0,'-black','LineWidth',1);    
for k=1:aSigCount
    aAmpl = aSource(:,k);  
    if aNorm == 1
        maxd = max(aAmpl);
    else
        maxd = 1;
    end;
    hold on;    
    plot(k+aAmax*(aAmpl/maxd),aTime,'Color',gwlGetColor(3),'LineStyle','-','LineWidth',0.5);     
    hold off;
end;
axis([0,aSigCount+1,aTmin,aTmax]);
set(gca,'LineWidth',0.5);
set(gca,'FontSize',aFont.Size-2);
set(gca,'Box','Off');
set(gca,'YDir','reverse');

if strcmp(aXLab,'')==1
    set(gca,'XTickLabel',{});
else
    gwlLabel('X',aXLab);
end;

if strcmp(aYLab,'')==1
    set(gca,'YTickLabel',{});
else
    gwlLabel('Y',aYLab);
end;

if (nargin > 12) 
    aSx = (aSigCount-1)/100;
    aSy = aTmin+(aTmax-aTmin)/100;
    gwlText(aSx,aSy,aText);
end;

if (nargin > 13) 
    gwlTitle(aTitle);
end;

grid on;
