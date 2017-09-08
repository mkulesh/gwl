function map = gwlhsv(m)
% map = gwlhsv(m)
% This is a modification of standart HSV color map. In this modification
% we use white color for minimal and maximal values of surface plot
%
% This is an internal procedure of GWL mshell library and is used to make colors 
% of all GWL plots uniform.
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

if nargin < 1, m = size(get(gcf,'colormap'),1); end
h = (0:m-1)'/max(m,1);
if isempty(h)
  map = [];
else
  map = hsv2rgb([h ones(m,2)]);
end
map(1,:) = [1 1 1];
map(m,:) = [1 1 1];

