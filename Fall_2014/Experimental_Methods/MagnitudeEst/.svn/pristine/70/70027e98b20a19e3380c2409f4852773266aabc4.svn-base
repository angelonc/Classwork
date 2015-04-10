function response = DoATrialME(win, params, lineInfo, standardInfo)
% DoATrialME - Run a trial of the magnitude estimation experiment.
%
% Syntax:
% response = DoATrialME(win, params, lineInfo, standardInfo)
%
% Input:
% win (GLWindow) - The GLWindow object.
% params (struct) - Structure containing the experimental parameters
%   defined in MagEstClass.m.
% lineInfo (struct) - Structure containing the line information.
% standardInfo (struct) - Structure containing the standard information.
%
% Output:
% response (scalar) - Trial response.
%
% NOTE FOR IMPROVEMENT: This routine could do a better job of checking
% for garbage input.

if nargin ~= 4
	error(help('DoATrialME'));
end

% If standard on every trial.
if params.standardAlways
	% Flush the mouse click buffer.
	mglListener('getAllMouseEvents');
	
	% Render the standard line.
	ShowTheLineME(win, params, standardInfo);
	
	if params.trialDuration == -1
		% Wait until the use clicks the mouse button.
		mglGetMouseEvent(Inf);
	else
		mglWaitSecs(params.standardDuration/1000);
	end
	
	% Turn the line off.
	win.disableAllObjects;
	win.draw;
	
	% Wait the standard ITI.
	mglWaitSecs(params.standardITI/1000);
end

% Show the trial
mglListener('getAllMouseEvents'); % Flush the mouse buffer.
ShowTheLineME(win, params, lineInfo);
if params.trialDuration == -1
	mglGetMouseEvent(Inf);
else
	mglWaitSecs(params.trialDuration/1000);
end
win.disableAllObjects;
win.draw;

% Add the response text prompts and wait for a response.
textTop = 250;
lineOffset = 60;
mglGetKeyEvent;
win.addText('"q" to quit"', 'Center', [0 textTop - lineOffset], 'Name', 'prompt1', 'FontSize', 80);
win.addText('"r" to repeat"', 'Center', [0 textTop - 2*lineOffset], 'Name', 'prompt2', 'FontSize', 80);
win.addText('"s" to show standard"', 'Center', [0 textTop - 3*lineOffset], 'Name', 'prompt3', 'FontSize', 80);
mglWaitSecs(params.blankDuration/1000);

response = NaN;
while isnan(response)
	response = mglGetEchoString(win, 'Enter magnitude estimate: ', [0 0], [1 1 1], 80);
	
	switch response
		case 'q'
			response = -1;
			
		case 's'
			response = -2;
			
		case 'r'
			response = -3;
			
		otherwise
			% Convert the string to a number.
			response = str2double(response);
	end
end
