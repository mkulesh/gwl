function gwlLabel(aLab, aText)
% gwlLabel(aLab, aText)
% This procedure puts a label aText to the x-axis (aLab='X'),
% y-axis (aLab='Y') or z-axis (aLab='Z') of a subplot.
%
% This is an internal procedure of GWL mshell library and is used in plotting
% methods like gwlPlotFunction, gwlPlotImage etc. to make labels of all 
% GWL plots uniform. 
%
% This file is part of the GWL library. Copyright (C) 2006-2007 Mikhail Kulesh, 
% mkulesh@math.uni-potsdam.de

aFont = gwlGetFont;
if aLab == 'X'
    xlabel(aText,'FontName',aFont.Name,'FontSize',aFont.Size,'FontWeight',aFont.Weight);
else if aLab == 'Y'
         ylabel(aText,'FontName',aFont.Name,'FontSize',aFont.Size,'FontWeight',aFont.Weight);
     else if aLab == 'Z'
             zlabel(aText,'FontName',aFont.Name,'FontSize',aFont.Size,'FontWeight',aFont.Weight);
          end;
     end;
end;
