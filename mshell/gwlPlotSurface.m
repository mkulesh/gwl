function gwlPlotSurface(aTime,aFreq,aSpec, aType,aPrec,aLineCount, aX1,aY1,aWidth,aHeight, aXLab,aYLab,aText,aTitle)
% gwlPlotImage(aTime,aFreq,aSpec, aX1,aY1,aWidth,aHeight, aXLab,aYLab [,aText,aTitle])
% This procedure summarize some standard procedures to plot a complex 2D array as a surface 
%
% Input parameters define the function and its changing interval, the window position of the plot as well as the labels:
%   aTime is an array corresponds to the horizontal axis
%   aFreq is an array corresponds to the vertical axis
%   aSpec is a 2-D array corresponds to the complex function values
%   aType defines the transformation rules:
%      aType=0: will be plotted abs(aSpec);
%      aType=1: will be plotted atan2(aSpec) using aPrec as a cut-off;
%      aType=2: will be plotted aSpec (only if aSpec is real)
%   aLineCount defines the number of lines for the surface (see contourf for details)
%   aX1,aY1,aWidth,aHeight describe the window position of the plot (see subplot for details)
%   aXLab, aYLab are a string variables correspond to the X- and Y- labels of the plot
%   aText is a string variable defines the text label in the left-top corner of the plot
%   aTitle is a string variable defines the plot title (see gwlTitle for details)
% 
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aFont = gwlGetFont;
aTmin = aTime(1);
aTmax = aTime(length(aTime));
aFmin = aFreq(1);
aFmax = aFreq(length(aFreq));

subplot('Position',[aX1,aY1,aWidth,aHeight])   
if aType == 0
    [C,h] = contourf(aTime,aFreq,abs(aSpec),aLineCount);
end;
if aType == 1
    [aSizey,aSizex]=size(aSpec);
    aMax = max(max(abs(aSpec)))*aPrec;
    for k=1:aSizex
        for m=1:aSizey
            if abs(aSpec(m,k)) > aMax
                aArg(m,k) = atan2(imag(aSpec(m,k)),real(aSpec(m,k)));
            else
                aArg(m,k) = -pi;
            end;
        end;
    end;
    [C,h] = contourf(aTime,aFreq,aArg,aLineCount);
end;
if aType == 2
    [C,h] = contourf(aTime,aFreq,aSpec,aLineCount);
end;
axis([aTmin,aTmax,aFmin,aFmax]);
set(gca,'LineWidth',1);
set(gca,'FontSize',aFont.Size-2);
set(gca,'Box','Off');

if(gwlGetColorShema == 0) colormap gwlgray; end;
if(gwlGetColorShema == 1) colormap gwlhot; end;
for k=1:length(h)
    set(h(k),'LineStyle','none');
end;

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
    aSy = aFmax-(aFmax-aFmin)/100;
    gwlText(aSx,aSy,aText);
end;

if (nargin > 13) 
    gwlTitle(aTitle);
end;

grid on;

