function aFont = gwlGetFont
% aFont = gwlGetFont
% This procedure returns font used in GWL
%
% This is an internal procedure of GWL mshell library and is used in plotting
% methods like gwlPlotFunction, gwlPlotImage etc. to make the font of all 
% GWL plots uniform.
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aFont.Size = 10;
aFont.Name = 'Times';
aFont.Weight = 'bold';

