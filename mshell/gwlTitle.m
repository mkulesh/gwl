function gwlTitle(aText)
% gwlTitle(aText)
% This procedure puts a title aText on the subplot panel. 
%
% This is an internal procedure of GWL mshell library and is used in plotting
% methods like gwlPlotFunction, gwlPlotImage etc. to make the titles of all 
% GWL plots uniform.
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aFont = gwlGetFont;
title(aText,'FontName',aFont.Name,'FontSize',aFont.Size,'FontWeight',aFont.Weight);
