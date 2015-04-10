% LineAdjust - Method of adjustment for line length comparison.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXPERIMENTAL PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.standardLength = 100;				% Standard line length (pixels)
params.standardArrowLength = 30;			% Standard arrow length (pixels)
params.standardArrowAngle = 135;			% Standard arrow angle (degrees re. inward)

params.stimLengths = [80 100 120];			% Starting lengths (pixels)
params.stimArrowLength = 30;				% Arrowhead lengths (pixels)
params.stimArrowAngle = 45;					% Arrowhead angles (degrees re. inward)

params.nBlocks = 1;							% Number of blocks
params.penSize = 4;							% Line thickness (pixels)
params.ITI = 100;							% Intertrial interval (millisecs) 
params.bgColor = [0 0 0];					% Color of background (RGB, each between 0 and 1)
params.lineColor = [0 1 0];					% Color of arrows (RGB, each between 0 and 1)
params.arrowColor = [1 1 1];				% Color of arrowheads (RGB, each between 0 and 1)
params.horizontalOffset = 100;				% Horizontal offset between standard and stimulus
params.jitterPosition = 100;				% Vary line position trial to trial?  0 means no.
											% Positive number specifies size of random jitter.
									
params.experimenter = 'dummy';              % Name of the experimenter.
params.subject = 'subject1';                % Name of the subject.
params.experimentName = 'LineAdjust';       % Root name of the experiment and data file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call driver
params.trialDuration = -1;					% This is ignored in this version.  Don't change.
LineAdjustDrvr(params);
