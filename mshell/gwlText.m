function gwlText(aX, aY, aText, aRotation)
% gwlText(aX, aY, aText [, aRotation])
% This procedure puts a label aText with coordinates aX and
% aY to the subplot panel. aRotation defines the direction of the text.
%
% This is an internal procedure of GWL mshell library and is used in plotting
% methods like gwlPlotFunction, gwlPlotImage etc. to make inscription of all 
% GWL plots uniform. 
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aFont = gwlGetFont;
if (nargin == 3) 
    aRotation = 0;
end;

text(aX,aY,aText,'FontName',aFont.Name,'FontSize',aFont.Size,'FontWeight',aFont.Weight,'VerticalAlignment','top','Rotation',aRotation);

