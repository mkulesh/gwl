function map = gwlgray(m)
% map = gwlgray(m)
% This procedure produces a modified gray-scaled color map that is used in gwlPlotImage
% and gwlPlotSurface procedures if gray-scaled palette is defined in GWL.
%
% This is an internal procedure of GWL mshell library and is used to make colors 
% of all GWL plots uniform.
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

if nargin < 1, m = size(get(gcf,'colormap'),1); end
g = (0:m-1)'/(max(m-1,1));
map = [1-g.*g 1-g.*g.*g 1-g.*g.*g];
