function HorizVertExp
% HorizVertExp - Program to collect magnitude estimates for line length.
%
% Description:
% Parameters here are set for measuring the horizontal
% vertical illusion using magnitude estimation.
%
% The lines are displayed both horizontally and vertically
% and estimates are collected separately for each.
%
% The observer's task is to enter the judged magnitude of
% each displayed line.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXPERIMENTAL PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.standardLength = 100;					  % Standard line length
params.stimLengths = [70 80 90 100 110 120 130];  % Test lengths (pixels)
params.stimAngles = [0 90];					      % Test angles (degrees re. horizontal)
params.arrowLengths = [20];                       % Arrowhead lengths (pixels)
params.arrowAngles = [0];						  % Arrowhead angles (degrees re. inward)
params.nBlocks = 1;								  % Number of blocks
params.penSize = 4;							      % Line thickness (pixels)
params.trialDuration = 500;					      % Trial duration (milliseconds) ...
												  % trialDuration = -1 means mouse click ends trial
params.blankDuration = 10;				          % Time before response box goes up (millisecs);
params.ITI = 100;								  % Intertrial interval (millisecs) 
params.bgColor = [0 0 0];						  % Color of background (RGB, each between 0 and 1)
params.lineColor = [1 1 1];		                  % Color of arrows (RGB, each between 0 and 1)
params.arrowColor = [1 1 1];		              % Color of arrowheads (RGB, each between 0 and 1)

params.jitterHPosition = 0; 					  % Vary horizontal line position trial to trial?
params.jitterVPosition = 0; 					  % Vary vertical line position trial to trial?
                                                  % 0 means no., nonzero gives
                                                  % maximum number of pixels of shift
                                                  % (left or right, up or down)
                                        
params.standardAlways = 0;						  % Show standard on every trial?
params.standardDuration = 750;				      % If standard shown, duration
params.standardITI = 200;						  % If standard shown, time to wait before line.

params.experimenter = 'dummy';                    % Experimenter
params.subject = 'subject1';                      % Name of the subject.
params.experimentName = 'HorizVertExp';		      % Root name of the experiment and data file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the actual experiment
MagnitudeEstDrvr(params);

% Gets rid of some Matlab syntax warnings.
%#ok<*NBRAK>
