function gwlPlotImage(aTime,aFreq,aSpec, aX1,aY1,aWidth,aHeight, aXLab,aYLab,aText,aTitle)
% gwlPlotImage(aTime,aFreq,aSpec, aX1,aY1,aWidth,aHeight, aXLab,aYLab [,aText,aTitle])
% This procedure summarize some standard procedures to plot a real 2D-array as an image 
%
% Input parameters define the function and its changing interval, the window position of the plot as well as the labels:
%   aTime is an array corresponds to the horizontal axis
%   aFreq is an array corresponds to the vertical axis
%   aSpec is a 2-D array corresponds to the real 2D-array
%   aX1,aY1,aWidth,aHeight describe the window position of the plot (see subplot for details)
%   aXLab, aYLab are a string variables correspond to the X- and Y- labels of the plot
%   aText is a string variable defines the text label in the left-top corner of the plot
%   aTitle is a string variable defines the plot title (see gwlTitle for details)
% 
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aFont = gwlGetFont;
aTmin = min(aTime);
aTmax = max(aTime);
aFmin = min(aFreq);
aFmax = max(aFreq);

subplot('Position',[aX1,aY1,aWidth,aHeight])   
imagesc(aTime,aFreq,aSpec);
axis([aTmin,aTmax,aFmin,aFmax]);
set(gca,'LineWidth',1);
set(gca,'FontSize',aFont.Size-2);
set(gca,'Box','Off');
set(gca,'YDir','normal');

if(gwlGetColorShema == 0) colormap gwlgray; end;
if(gwlGetColorShema == 1) colormap gwlhot; end;

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

if (nargin > 9) 
    aSx = aTmin+(aTmax-aTmin)/100;
    aSy = aFmax-(aFmax-aFmin)/100;
    gwlText(aSx,aSy,aText);
end;

if (nargin > 10) 
    gwlTitle(aTitle);
end;

grid on;
