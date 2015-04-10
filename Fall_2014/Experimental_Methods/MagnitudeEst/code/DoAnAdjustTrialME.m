function response = DoAnAdjustTrialME(win, params, lineInfo, standardInfo)
% DoAnAdjustTrialME - Run a trial of line adjustment.
%
% Syntax:
% response = DoAnAdjustTrialME(win, params, lineInfo, standardInfo)
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

% Show the trial.
ShowTheLineME(win, params, [standardInfo lineInfo])

% Loop until we get a valid response.
response = [];
lineDelta = 0;
adjustDelta = 5;
mglGetKeyEvent;
while isempty(response);
	key = mglGetKeyEvent(Inf);
	
	switch lower(key.charCode)
		% Abort.
		case 'q'
			response = -1;
			
		% Accept adjustment.
		case char(13)
			response = lineDelta;
			
		% Down arrow pressed.
		case char(31)
			lineDelta = lineDelta + adjustDelta;
			lineInfo.to(2) = lineInfo.to(2) + adjustDelta;
			lineInfo.upArrowTo(2) = lineInfo.upArrowTo(2) + adjustDelta;
			lineInfo.downArrowTo(2) = lineInfo.downArrowTo(2) + adjustDelta;
			ShowTheLineME(win, params, [standardInfo lineInfo]);
			
		% Up arrow pressed.
		case char(30)
			lineDelta = lineDelta - adjustDelta;
			lineInfo.to(2) = lineInfo.to(2) - adjustDelta;
			lineInfo.upArrowTo(2) = lineInfo.upArrowTo(2) - adjustDelta;
			lineInfo.downArrowTo(2) = lineInfo.downArrowTo(2) - adjustDelta;
			ShowTheLineME(win, params, [standardInfo lineInfo]);
			
		% Right arrow pressed.  Bigger adjustment size.
		case char(29)
			adjustDelta = 5;
			
		% Left arrow pressed.  Smaller adjustment size.
		case char(28)
			adjustDelta = 1;
	end	
end

% Blank the screen.
win.disableAllObjects;
win.draw;
