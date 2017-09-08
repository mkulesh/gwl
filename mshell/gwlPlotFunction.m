function gwlPlotFunction(aTime, aSource, aX1,aY1,aWidth,aHeight, aTmin,aTmax,aAmin,aAmax, aXLab,aYLab,aText,aTitle, aLog)
% gwlPlotFunction(aTime, aSource, aX1,aY1,aWidth,aHeight, aTmin,aTmax,aAmin,aAmax, aXLab,aYLab [,aText,aTitle,aLog])
% This procedure summarize some standard procedures to plot a function 
%
% Input parameters define the function and its changing interval, the window position of the plot as well as the labels:
%   aTime is an array corresponds to the function argument
%   aSource is an array corresponds to the function values
%   aX1,aY1,aWidth,aHeight describe the window position of the plot (see subplot for details)
%   aTmin,aTmax,aAmin,aAmax define the interval of function and its argument (see set(gca,'Xlim') for details)
%   aXLab, aYLab are a string variables correspond to the X- and Y- labels of the plot
%   aText is a string variable defines the text label in the left-top corner of the plot
%   aTitle is a string variable defines the plot title (see gwlTitle for details)
%   aLog is a string variable; if aLog == 'logx' then the logarithmic X-axes will be used,
%      if aLog == 'logy' then the logarithmic Y-axes will be used (see semilogx and semilogy for details)
% 
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de
             
aFont = gwlGetFont;

subplot('Position',[aX1,aY1,aWidth,aHeight]);
if(nargin > 14)
    if aLog == 'logx'
        semilogx(aTime,aSource,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    else
        semilogy(aTime,aSource,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
    end;
else
    plot(aTime,aSource,'Color',gwlGetColor(0),'LineStyle','-','LineWidth',1);
end;
set(gca,'Xlim',[aTmin aTmax]);
set(gca,'Ylim',[aAmin,aAmax]);
set(gca,'LineWidth',0.5);
set(gca,'FontSize',aFont.Size-2);
set(gca,'Box','Off');

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
    aSx = aTmin+(aTmax-aTmin)/100;
    aSy = aAmax-(aAmax-aAmin)/100;
    gwlText(aSx,aSy,aText);
end;

if (nargin > 13) 
    gwlTitle(aTitle);
end;

grid on;

