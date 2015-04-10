function LinePSE
% LinePSE
%
% Program to measure psychometric function for line length
% comparison.
%
% 5/25/98   dhb   Wrote it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXPERIMENTAL PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.standardLength = 100;					% Standard line length (pixels)
params.standardArrowLength = 30;			    % Standard arrow length (pixels)
params.standardArrowAngle = 135; 			    % Standard arrow angle (degrees re. inward)

params.stimLengths = [100 105 110 115 120 125 130 135 140]; % Test lengths (pixels)
params.stimArrowLength = 30;					% Arrowhead lengths (pixels)
params.stimArrowAngle = 45;					    % Arrowhead angles (degrees re. inward)

params.nBlocks = 8;								% Number of blocks
params.penSize = 4;								% Line thickness (pixels)
params.trialDuration = -1;						% Trial duration (milliseconds) ...	
												% -1 means wait for response
params.ITI = 100;								% Intertrial interval (millisecs) 
params.bgColor = [0 0 0];						% Color of background (RGB, each between 0 and 1)
params.lineColor = [0 1 0];				        % Color of arrows (RGB, each between 0 and 1)
params.arrowColor = [1 1 1];		            % Color of arrowheads (RGB, each between 0 and 1)
params.horizontalOffset = 100;				    % Horizontal offset between standard and stimulus
params.jitterPosition = 100; 				    % Vary line position trial to trial?  0 means no.
											    % Positive number specifies size of random jitter.

params.experimenter = 'dummy';                  % Name of the experimenter.
params.subject = 'subject1';                    % Name of the subject.
params.experimentName = 'LinePSE';              % Root name of the experiment and data file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the driver.
LinePSEDrvr(params);
