function AddMagEstPath
% AddMagEstPath - Dynamically adds the MagnitudeEst folder to the path.
%
% Syntax:
% AddMagEstPath

p = fileparts(fileparts(which('AddMagEstPath')));
if isempty(strfind(path, p))
	addpath(RemoveSVNPaths(genpath(p)), '-end');
end
