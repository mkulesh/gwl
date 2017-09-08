function gwlPlotPngImage(aImage, aX1,aY1,aWidth,aHeight, aText,aTitle)
% gwlPlotPngImage(aImage, aX1,aY1,aWidth,aHeight [, aText,aTitle])
% This procedure summarize some standard procedures to plot an image from a file 
%
% Input parameters define the function and its changing interval, the window position of the plot as well as the labels:
%   aImage is the name of the graphical file, for example, in PNG-format
%   aX1,aY1,aWidth,aHeight describe the window position of the plot (see subplot for details)
%   aText is a string variable defines the text label in the left-top corner of the plot
%   aTitle is a string variable defines the plot title (see gwlTitle for details)
% 
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

[aImage] = imread(aImage);
subplot('Position',[aX1,aY1,aWidth,aHeight]);
image(aImage);

set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});
set(gca,'Box','Off');
set(gca,'Visible','Off');

if (nargin > 5) 
    aTmin = 1;
    aAmax = 1;
    [aTmax,aAmin,aCol] = size(aImage);
    aSx = aTmin+(aTmax-aTmin)/100;
    aSy = aAmax-(aAmax-aAmin)/100;
    gwlText(aSx,aSy,aText);
end;

if (nargin > 6) 
    gwlTitle(aTitle);
end;
