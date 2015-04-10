function LineAdjustDrvr(params)
% LineAdjustDrvr - Driver program for line PSE experiment.
%
% Syntax:
% LineAdjustDrvr(params)

% Make sure we're using the right version of GLWindow.
UseClassesDev;

% Convenience parameters
nStimLengths = length(params.stimLengths);
nIndVarCols = 1;
nDataRows = nStimLengths;
nDataColumns = nIndVarCols + params.nBlocks;
nTrials = params.nBlocks*nDataRows;

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
%		Column 1:																Stimulus lengths
%		Columns 5:(nBlocks+4)										Estimates
%
% On abort, any unfilled entries are still at -1.
theData = -1*ones(nDataRows, nDataColumns);
index = 1;
for i = 1:nStimLengths
	theData(index,1) = params.stimLengths(i);
	index = index+1;
end

% Seed the random number generator
ClockRandSeed;

% Open up the experimental window and set globals.
win = OpenEXPWindowME(params.bgColor);

try
	ListenChar(2);
	
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
	% combinations and display.
	for i = 1:params.nBlocks
		% Shuffle trials
		indices = Shuffle(1:nDataRows);
		
		% Do the trials
		for j = 1:nDataRows
			index = indices(j);
			dataRow = index;
			dataColumn = nIndVarCols+i;
			origLength = theData(dataRow,1);
			
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
			
			response = DoAnAdjustTrialME(win, params, tempInfos, standardInfo);
			theLength = origLength + response;
			
			% This will cause execution to abort the program immediately.
			if response == -1
				error('abort');
			end
			
			% Save the data
			theData(dataRow, dataColumn) = theLength;
			
			% Wait intertrial interval
			mglWaitSecs(params.ITI/1000);
		end
	end
	
	% Get total time
	totalSecs = mglGetSecs - startRunSecs;
	
	% Figure out some data saving parameters.
	dataFolder = sprintf('%s/data/%s/%s/%s', fileparts(fileparts(which('MagEstClass'))), ...
		params.experimenter, params.experimentName, params.subject);
	if ~exist(dataFolder, 'dir')
		mkdir(dataFolder);
	end
	dataFile = sprintf('%s/%s-%d.csv', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.csv'));

	% Write out CSV data.
	c = CSVFile(dataFile, true);
	dataColumns = {'Starting Length'};
	for i = 1:params.nBlocks
		dataColumns{end+1} = sprintf('Response %d', i); %#ok<AGROW>
	end
	for i = 1:length(dataColumns)
		c = c.addColumn(dataColumns{i}, 'd');
		c = c.setColumnData(dataColumns{i}, theData(:,i));
	end
	c.write;
	
	% Print out mean adjustment.
	WaitSecs(0.5);
	allAdjust = theData(:, 2:1 + params.nBlocks);
	allAdjust = allAdjust(:);
	nAdjust = length(allAdjust);
	allMean = mean(allAdjust);
	allSEM = std(allAdjust) / sqrt(nAdjust);
	fprintf('Mean adjustment = %g +/- %g\n', allMean, allSEM);
	matchString = sprintf('Mean adjustment = %g +/- %g\n', allMean, allSEM);
	
	% Little screen display of results.
	mglGetKeyEvent;
	matchInfo = PrecomputeLineME(allMean, 90, params.stimArrowLength, ...
		params.stimArrowAngle, params.horizontalOffset/2);
	standardInfoN = PrecomputeLineME(params.standardLength, 90, ...
		params.standardArrowLength, 0, -params.horizontalOffset/2);
	matchInfoN = PrecomputeLineME(allMean, 90, params.stimArrowLength, 0,...
		params.horizontalOffset/2);
	ShowResultLine(win, params, [standardInfo matchInfo], matchString);

	while true
		key = mglGetKeyEvent(Inf);
		
		switch lower(key.charCode)
			case 'a'
				ShowResultLine(win, params, [standardInfo matchInfo], matchString);
				
			case 'n'				
				ShowResultLine(win, params, [standardInfoN matchInfoN], matchString);
				
			case 'q'
				break;
		end
	end
	
	% Close the experimental window
	ListenChar(0);
	win.close;
catch e
	ListenChar(0);
	win.close;
	
	if strcmp(e.message, 'abort')
		fprintf('User quit experiment early\n');
	else
		rethrow(e);
	end
end


function ShowResultLine(win, params, lineInfo, matchString)
% Make sure we don't render stuff we don't want.
win.disableAllObjects;

win.addText(matchString, 'Name', 't0', 'Center', [0 300], 'FontSize', 80);
win.addText('No Arrows - N', 'Name', 't1', 'Center', [0 250], 'FontSize', 80);
win.addText('With Arrows - A', 'Name', 't2', 'Center', [0 200], 'FontSize', 80);
win.addText('Quit - Q', 'Name', 't3', 'Center', [0 150], 'FontSize', 80);
						
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
