function behOrder(sub)

% This function generates the presentation order for a subjects
% test phase of the experiment
%
% What we want:
% - 6 blocks
% - all experimental objects + 16 foils (80 total)
% - each block shows different views (random)

if ~exist('sub','var') sub = 99; end

% Input
Dir = ['~/Documents/Classwork/Epstein_Rotation/LocalGeometry/' ...
             'Behavior/'];
fn = [Dir sprintf('Configs/sub%02d_config.mat',sub)];

% Output
outDir = [Dir sprintf('Data/%02d/',sub)];
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% Load object indexes
load(fn);

% Make a list of our test objects
foils   = index.names(index.conditions(:,1) == 1);
exp     = index.names(index.conditions(:,2) == 1);
objList = [foils; exp];

% Randomize presentation angle for all trials and blocks
for i = 1:length(objList)
    angle1(i,:) = shuffle([1 2 3 5 6 7]);
end

% Index object names for textures
for i = 1:length(objList)
    str = strread(objList{i},'%s','delimiter','_*.');
    if strcmp(str{1}(1),'f')
        f = 1;
    else
        f = 2;
    end
        
    objIdx(i,:) = [f str2num(str{1}(2:3))];
end

% Create randomized block order
for i = 1:size(angle1,2)
    % Shuffle idx
    idx = randperm(size(angle1,1))';
    
    % Shuffle the object list and angle list
    angle(:,i) = angle1(idx,i);
    order{:,i} = [objIdx(idx,:) angle(:,i)];
    names(:,i) = objList(idx);
end

beh.names = names;
beh.order = order;
beh.angle = angle;

% Write to file
save([outDir 'test_order.mat'],'beh');






