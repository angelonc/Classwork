function ShowStandardME(win, params, standardInfo)
% ShowStandardME(win, params, standardInfo)
%
% Description:
% Show the standard and wait for a response.
%
% Input:
% win (GLWindow) - The GLWindow object.
% params (struct) - Structure containing the experimental parameters
%   defined in MagEstClass.m.
% stanardInfo (struct) - Struct containing the line parameters.

if nargin ~= 3
	error('Usage: ShowStandardME(win, params, standardInfo)');
end

% Make sure we don't have anything else drawable in the queue.
win.disableAllObjects;

% Create the names for the text objects we'll use.  We'll delete them at
% then end of the function.
for i = 1:2
	textNames{i} = GenerateRandomString(20); %#ok<AGROW>
end

% Show the line and message.
win.addText('The standard below is 100.', 'Name', textNames{1}, 'Center', [0 50], 'FontSize', 80);
win.addText('Hit any key to continue.', 'Name', textNames{2}, 'Center', [0 100], 'FontSize', 80);
win.addLine(standardInfo.from, standardInfo.to, params.penSize, params.lineColor, 'Name', 'standard');
win.draw;

FlushEvents;
GetChar;

% Turn off the stuff on the screen.
for i = 1:length(textNames)
	win.deleteObject(textNames{i});
end
win.disableAllObjects;
win.draw;
