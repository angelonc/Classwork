function LinePSEDrvr(params)
% LinePSEDrvr - Driver program for line PSE experiment.
%
% Syntax:
% LinePSEDrvr(params)

% Use the correct version of Classes
UseClassesDev;

% Make sure 'params' was passed.
if nargin ~= 1
	error(help('LinePSEDrvr'));
end

% Seed the random number generator
ClockRandSeed;

ListenChar(2);

% Convenience parameters
nStimLengths = length(params.stimLengths);
nIndVarCols = 1;
nDataRows = nStimLengths;
nDataColumns = nIndVarCols + params.nBlocks + 2;
nTrials = params.nBlocks * nDataRows;

% % Load in calibration infomation if it isn't there
% DO_CALIBRATED = 0;
% if (DO_CALIBRATED)
% 	if (exist('calME') ~= 1)
% 		disp('Loading calibration file');
% 		calME = LoadCalFile(whichScreen);
% 		calME = SetGammaMethod(calME,1);
% 	else
% 		disp('Using calibration file in memory');
% 	end
% 	bgLin = DeviceToSettings(calME,bgColor'/255)';
% 	lineLin = DeviceToSettings(calME,lineColor'/255)';
% 	arrowLin = DeviceToSettings(calME,arrowColor'/255)';
% else
% 	bgLin = bgColor;
% 	lineLin = lineColor;
% 	arrowLin = arrowColor;
% end

% Create the data array.  Format is:
%		Column 1:				Stimulus lengths
%		Columns 5:(nBlocks+4)	Estimates
%
% On abort, any unfilled entries are still at -1.
theData = -1*ones(nDataRows,nDataColumns);
index = 1;
for i = 1:nStimLengths
	theData(index,1) = params.stimLengths(i);
	index = index+1;
end

% Open the experimental window.
win = OpenEXPWindowME(params.bgColor);

try
	% Precompute line display points
	for i = 1:nDataRows
		lineInfos(i) = PrecomputeLineME(theData(i,1), 90, params.stimArrowLength, ...
			params.stimArrowAngle, params.horizontalOffset/2); %#ok<AGROW>
	end
	standardInfo = PrecomputeLineME(params.standardLength, 90, params.standardArrowLength, ...
		params.standardArrowAngle, -params.horizontalOffset/2);
	
	% Get start of run
	startRunSecs = mglGetSecs;
	
	% Do the trials.  In each block, we shuffle all possible
	% combinations and display
	for i = 1:params.nBlocks
		% Shuffle trials
		indices = Shuffle(1:nDataRows);
		
		% Do the trials
		for j = 1:nDataRows
			index = indices(j);
			dataRow = index;
			dataColumn = nIndVarCols + i;
			
			% Do a single trial and record response.  Abort on response of -1.
			tempInfos = lineInfos(index);
			if (params.jitterPosition > 0)
				theJitter = Ranint(2*params.jitterPosition) - params.jitterPosition;
				
				% Add the jitter to all the horizontal components of the line
				% info.
				fieldNames = fieldnames(tempInfos);
				for k = 1:length(fieldNames)
					tempInfos.(fieldNames{k})(1) = tempInfos.(fieldNames{k})(1) + theJitter;
				end
			end
			
			response = DoAPSETrialME(win, params, tempInfos, standardInfo);
			if response == -1
				error('abort');
			end
			
			% Save the data
			theData(dataRow,dataColumn) = response;
			
			% Wait intertrial interval
			mglWaitSecs(params.ITI/1000);
		end
	end
	
	% Get total time
	totalSecs = mglGetSecs - startRunSecs;
	
	% Close the experimental window
	win.close;
	
	% Analyze data if there is enough
	% Fit a logistic through the data.
	pRight = sum(theData(:, nIndVarCols+1:nIndVarCols+params.nBlocks)')' / params.nBlocks;
	[a, b] = FitLogistic(theData(:,1), pRight);
	fitLogistic = ComputeLogistic(theData(:,1), a, b);
	p10 = InvertLogistic(0.10, a, b);
	p25 = InvertLogistic(0.25, a, b);
	p50 = InvertLogistic(0.50, a, b);
	p75 = InvertLogistic(0.75, a, b);
	p90 = InvertLogistic(0.90, a, b);
	thresh75 = (p75-p25)/2;
	thresh90 = (p90-p10)/2;
	
	% Figure out some data saving parameters.
	ci = GetComputerInfo;
	experimenter = ci.userShortName;
	dataFolder = sprintf('%s/data/%s/%s/%s', fileparts(fileparts(which('MagEstClass'))), ...
		experimenter, params.experimentName, params.subject);
	if ~exist(dataFolder, 'dir')
		mkdir(dataFolder);
	end
	dataFile = sprintf('%s/%s-%d.csv', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.csv'));
	
	% Write the CSV file.
	c = CSVFile(dataFile, true);
	dataColumns = {'Line Length'};
	for i = 1:params.nBlocks
		dataColumns{end+1} = sprintf('Response %d', i); %#ok<AGROW>
	end
	dataColumns(end+1:end+2) = {'pRight', 'Logistic Fit'};
	theData = [theData, pRight, fitLogistic];
	for i = 1:length(dataColumns)
		c = c.addColumn(dataColumns{i}, 'g');
		c = c.setColumnData(dataColumns{i}, theData(:,i));
	end
	c.write;
	
	% Put up a plot of the data
	hold off;
	plot(theData(:,1), pRight, 'r+');
	hold on;
	plot(theData(:,1), fitLogistic, 'g');
	axis([params.stimLengths(1) params.stimLengths(nDataRows) 0 1]);
	if exist('THRESHOLD', 'var')
		if THRESHOLD
			title(sprintf('75%% threshold = %g', p75));
		else
			title(sprintf('PSE = %g', p50));
		end
	else
		title(sprintf('PSE = %g', p50));
	end
	hold off;
	
	% Print the PSE to window as well
	fprintf('PSE = %g\n', p50);
	fprintf('75 percent threshold = %g\n', p75);
	fprintf('90 percent threshold = %g\n', p90);
	fprintf('Logistic parameters: a = %g, b = %g\n', a, b);
	fprintf('Run took %g minutes\n\n', totalSecs / 60);
	
	% Close the experimental window and turn off character listening.
	ListenChar(0);
	win.close;
catch e
	ListenChar(0);
	win.close;
	
	if strcmp(e.message, 'abort')
		fprintf('Experiment aborted by user.\n');
	else
		rethrow(e);
	end
end




