% MagEstExampleAnalysis - Example analysis script for MagEst programs.

% Clear out old stuff
clear; close all;

% Define parameters that describe the experiment that was run.
experimenter = 'dummy';
experimentName = 'MullerLyer';
subject = 'subject1';

% List of file numbers we want to analyze.  These are the numbers appended
% to each of the data files.
fileNumbers = [1 2];

% Now we can extract the data into a MxNxP matrix along with the headers
% for each column of data.  Each row represents a particular stimulus and
% all its measurements.  Each column represents a particular
% parameter/measurement for a set of stimuli.  P represents a set of data
% for a given run of the experiment or data file.
[data, columnHeaders] = MagEstExtractData(experimenter, experimentName, subject, fileNumbers);
