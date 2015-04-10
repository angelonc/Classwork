function ShowTheLineME(win, params, lineInfo)
% ShowTheLineME - Draw line stimuli.
%
% Syntax:
% ShowTheLineME(win, params, lineInfo)
%
% Description:
% Draws line stimuli about the center of the display.  All objects added
% to the GLWindow prior to this call will be disabled.
%
% Input:
% win (GLWindow) - The GLWindow object.
% params (struct) - Structure containing the experimental parameters
%   defined in MagEstClass.m.
% lineInfo (1xN struct) - Struct array containing line information to
%   render.

if nargin ~= 3
	error(help('ShowTheLineME'));
end

% Make sure we don't render stuff we don't want.
win.disableAllObjects;
						
% Draw the line
n = length(fieldnames(lineInfo));
for i = 1:length(lineInfo)
	win.addLine(lineInfo(i).from, lineInfo(i).to, params.penSize, params.lineColor, 'Name', sprintf('l_%d', (i-1)*n + i + 4));
	win.addLine(lineInfo(i).from, lineInfo(i).upArrowFrom, params.penSize, params.arrowColor, 'Name', sprintf('l_%d', (i-1)*n + i));
	win.addLine(lineInfo(i).from, lineInfo(i).downArrowFrom, params.penSize, params.arrowColor, 'Name', sprintf('l_%d', (i-1)*n + i + 1));
	win.addLine(lineInfo(i).to, lineInfo(i).upArrowTo, params.penSize, params.arrowColor, 'Name', sprintf('l_%d', (i-1)*n + i + 2));
	win.addLine(lineInfo(i).to, lineInfo(i).downArrowTo, params.penSize, params.arrowColor, 'Name', sprintf('l_%d', (i-1)*n + i + 3));
end
win.draw;
