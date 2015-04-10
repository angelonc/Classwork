function response = DoAPSETrialME(win, params, lineInfo, standardInfo)
% DoAPSETrialME - Runs a trial of the line PSE experiment.
%
% Syntax:
% response = DoAPSETrialME(win, params, lineInfo, standardInfo)
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


% Show the trial
ShowTheLineME(win, params, [standardInfo lineInfo]);
if params.trialDuration > 0
	mglWaitSecs(params.trialDuration/1000);
	
	% Turn the lines off.
	win.disableAllObjects;
	win.draw;
end

% Loop until a valid response is given.
response = [];
FlushEvents;
while isempty(response)
	key = GetChar;
	
	switch lower(key)
		case 'q'
			response = -1;
			
		case 'd'
			response = 0;
			
		case 'k'
			response = 1;
	end
end

if params.trialDuration <= 0
	% Turn the lines off.
	win.disableAllObjects;
	win.draw;
end
