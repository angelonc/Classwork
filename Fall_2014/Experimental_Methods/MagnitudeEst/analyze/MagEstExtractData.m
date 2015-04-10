function [data, columnHeaders] = MagEstExtractData(experimenter, experimentName, subjectName, fileNumbers)
% MagEstExtractData - Extracts MagEst data from .csv files.
%
% Syntax:
% [data, columnHeaders] = MagEstExtractData(experimenter, experimentName, subjectName, fileNumbers)
%
% Description:
% Magnitude estimation data are saved as .csv files.  This function reads
% in a set of .csv files specified by the function arguments and returns a
% matrix containing the data along with the column headers for the data.
%
% Input:
% experimenter (string) - Experimenter's name as found in the data folder.
% experimentName (string) - Name of the experiment as found in the data
%   folder.
% subjectName (string) - Name of the subject as found in the data folder.
% fileNumbers (vector) - Vector of integer values that represent the data
%   files.  For instance, there may be datafile-1.csv, datafile-2.csv, etc.
%   To specify the first two files, you would pass the vector [1 2].  Any
%   other files ares specified in the same manner.
%
% Output:
% data (MxNxP) - Each row represents a stimulus.  Each column
%   represents a stimulus property or a trial result.  P represents a
%   particular run of the experiment.  P will be the same as the length of
%   the file numbers specified.
% columnHeaders (cell array) - The column headers for every column (N) of
%   the data matrix.

% Check the number of input arguments.
error(nargchk(4, 4, nargin));

if nargout ~= 2
	error('Make sure you are storing both return values of MagEstExtractData.');
end

% Construct the data file folder and make sure it exists.
dataDir = sprintf('%s/data/%s/%s/%s', fileparts(fileparts(which(mfilename))), experimenter, ...
	experimentName, subjectName);
if ~exist(dataDir, 'dir')
	error('Cannot find folder %s\n', dataDir); %#ok<SPERR>
end

% Make sure the fileNumbers vector is ok.
if ~isvector(fileNumbers)
	error('"fileNumbers" must be a vector.');
end
if any(fileNumbers < 1)
	error('"fileNumbers" must contain only positive integers.');
end

% Get a list of all the .csv files in the data folder.
d = dir(sprintf('%s/*.csv', dataDir));

% Filter out only the files specified.
dataFileNames = {};
matcher = sprintf('%d|', fileNumbers);
matcher = matcher(1:end-1);
pattern = sprintf('.*-(%s)\\.csv$', matcher);
for i = 1:length(d)
	if regexp(d(i).name, pattern)
		dataFileNames{end+1} = d(i).name; %#ok<AGROW>
	end
end

% The number of data files returned should match the number of file numbers
% specified.
if length(dataFileNames) ~= length(fileNumbers)
	s = sprintf('%d ', fileNumbers);
	error('Cannot find all the files specified by file numbers %s.\n', s(1:end-1));
end

% Open the first data file to see how big it is so we can define our data
% dimensions.
csvFile = CSVFile(sprintf('%s/%s', dataDir, dataFileNames{1}));
columnHeaders = csvFile.getColumnHeaders;
cData = csvFile.getColumnData(columnHeaders{1});
numRows = length(cData);
numCols = length(columnHeaders);

% Pre-allocate our data matrix.
data = zeros(numRows, numCols, length(dataFileNames));

% Extract each file's data and stick it into the data matrix.
for P = 1:length(dataFileNames)
	% Open the data file.
	csvFile = CSVFile(sprintf('%s/%s', dataDir, dataFileNames{P}));
	
	% Make sure it has the proper number of rows and columns.
	if (length(csvFile.getColumnHeaders) ~= numCols) || ...
	   (length(csvFile.getColumnData(columnHeaders{1})) ~= numRows)
		error('Data in data file %s doesn''t match the size of the other file data.');
	end
	
		for N = 1:numCols
			columnData = csvFile.getColumnData(columnHeaders{N});
			
			for M = 1:numRows
				data(M,N,P) = str2double(columnData{M});
			end
		end
	end
end
