function map = gwlhot(m)
% map = gwlhot(m)
% This procedure produces a white-blue-black color map that is used in gwlPlotImage
% and gwlPlotSurface procedures if gray-scaled palette is defined in GWL.
%
% This is an internal procedure of GWL mshell library and is used to make colors 
% of all GWL plots uniform.
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = fix(3/8*m);
r = [(1:n)'/n; ones(m-n,1)];
g = [zeros(n,1); (1:n)'/n; ones(m-2*n,1)];
b = [zeros(2*n,1); (1:m-2*n)'/(m-2*n)];
map = [1-r 1-g 1-b];
