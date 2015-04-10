function mazeOrder(sub)

dbstop if error

% This function makes the order of maze objects and indexes the
% objects for display in the test phase. It's really a shitshow.

if nargin < 1 sub = 99; end

% Attend to furniture for even sub #s, vehicles for odd
attend = [0 0];
attend(mod(sub,2) + 1) = 1;

% Offset values
posOff = .65;
angOff = 45;


%% I/O
% Directories
behDir = '~/Documents/Classwork/Epstein_Rotation/LocalGeometry/Behavior/';
unityDir = '~/Documents/Classwork/Epstein_Rotation/LocalGeometry/Maze_Unity/Assets/Config_Files/';

% Load base transforms
filePath = [behDir 'Configs/objectTransforms_base.txt'];
transforms = dlmread(filePath,',');
% plot(transforms(:,1),transforms(:,3))


%% Balance experimental conditions
DP     = [0 1];
junc   = [3 4];
turn   = [-1 1];
attend = [0 1];
pos    = [-1 1];

cnt = 1;
for i = 1:length(DP)
    for j = 1:length(junc)
        for k = 1:length(turn)
            for l = 1:length(attend)
                for m = 1:length(pos)
                    design(cnt,1) = DP(i);
                    if DP(i)
                        design(cnt,2) = junc(j);
                    else
                        design(cnt,2) = 0;
                    end
                    design(cnt,3) = turn(k);
                    design(cnt,4) = attend(l);
                    design(cnt,5) = pos(m);
                    cnt = cnt + 1;
                end
            end
        end
    end
end


%% Randomize and match to experimental spatial order for each maze
% Index objects (col 1 = furniture, col 2 = vehicles)
objNums = [(1:48)' (1:48)'];
tmp = [zeros(1,24) ones(1,24)];
% Assign objects
maze2 = [shuffle(tmp)' shuffle(tmp)'];


% Set patterns
dpPatt = [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]';
juncPatt(:,1) = [0 4 0 3 0 3 0 4 0 4 0 3 0 3 0 4 ...
                 4 0 3 0 3 0 4 0 4 0 3 0 3 0 4 0]';
juncPatt(:,2) = [0 3 0 4 0 4 0 3 0 3 0 4 0 4 0 3 ...
                 3 0 4 0 4 0 3 0 3 0 4 0 4 0 3 0]';
turnPatt = [repmat([-1 -1 1 1],1,4) repmat([1 1 -1 -1],1,4)]';

% For each maze
for i = 1:2
    %% Experimental Order
    order = nan(32,5);
    order(:,1:3) = [dpPatt juncPatt(:,i) turnPatt];
    
    % Randomize order
    idx = randperm(length(design));
    tmp = design(idx,:);
    
    % Get out the patterns
    patts = unique(order(:,1:3),'rows');
    for j = 1:length(patts)
        % Match patterns
        orderMatchIdx = find(sum(order(:,1:3) == repmat(patts(j,:),length(tmp),1),2)==3);
        designMatchIdx = find(sum(tmp(:,1:3) == repmat(patts(j,:),length(tmp),1),2)==3);
        
        % Assign to match
        order(orderMatchIdx,:) = tmp(designMatchIdx,:);
    end
    
        
    %% Object Locations and Transforms
    % Label experimental and non-experimental locations and overall turn index
    maze(i).idx.experimental  = [0 ones(1,16) 0 0 0 0 0 0 ones(1,16) 0]';
    maze(i).idx.turns = [0 repmat([-1 -1 1 1],1,4) ...
                        -1 -1 0 0 -1 -1 ...
                        repmat([1 1 -1 -1],1,4) 0]';
    
    % Randomize non-experimental object position offsets/attended
    nexpOff = shuffle([1 1 1 1 -1 -1 -1 -1]);
    nexpAttend = shuffle([1 1 1 1 0 0 0 0]);
    
    % Make matrix (each row = object location)
    objOrder = [NaN NaN 0 nexpAttend(1) nexpOff(1)];
    objOrder = [objOrder; order(1:length(order)/2,:)];
    objOrder = [objOrder; [nan(6,1) nan(6,1) [-1 -1 0 0 -1 -1]' ...
                        nexpAttend(2:7)' nexpOff(2:7)']];
    objOrder = [objOrder; order((length(order)/2) + 1:end,:)];
    objOrder = [objOrder; [NaN NaN 0 nexpAttend(end) nexpOff(end)]];
    
    
    % Add in offsets (sooo ugly)
    for j = 1:length(objOrder)
        % For left positioned objects
        if (objOrder(j,5) < 0)
            % For angles
            switch (transforms(j,end))
              case -90
                posOffs(j,:) = [0 0 -posOff];
                angOffs(j) = -angOff;
              case 90
                posOffs(j,:) = [0 0 posOff];
                angOffs(j) = -angOff;
              case 180
                posOffs(j,:) = [posOff 0 0];
                angOffs(j) = -angOff;
              case 0
                posOffs(j,:) = [-posOff 0 0];
                angOffs(j) = -angOff;
            end
        % For right positioned objects
        else
            % For angles
            switch (transforms(j,end))
              case -90
                posOffs(j,:) = [0 0 posOff];
                angOffs(j) = angOff;
              case 90
                posOffs(j,:) = [0 0 -posOff];
                angOffs(j) = angOff;
              case 180
                posOffs(j,:) = [-posOff 0 0];
                angOffs(j) = angOff;
              case 0
                posOffs(j,:) = [posOff 0 0];
                angOffs(j) = angOff;
            end
        end
    end
    
    % Calculate new transforms
    maze(i).transforms.original = transforms;
    maze(i).transforms.new = transforms + [posOffs angOffs'];

    % Concatenate
    objOrder = [objOrder maze(i).transforms.new];
    
    
    %% Object Order
    % Foils
    objs = reshape(objNums(maze2==i-1),24,2);
    foilIdx = [shuffle([zeros(20,1); ones(4,1)]) shuffle([zeros(20,1); ones(4,1)])];
    foils = reshape(objs(foilIdx == 1),4,2);
    cnt = 1;
    for j = 1:size(foils,2)
        for k = 1:size(foils,1)
            if j == 1
                str = 'f';
            else
                str = 'v';
            end
            foilObjects{cnt,1} = sprintf('%s%02d',str,foils(k,j));
            cnt = cnt + 1;
        end
    end
    
    % Non foils
    nFoils = reshape(objs(foilIdx == 0),20,2);
    for j = 1:size(nFoils,2)
        for k = 1:size(nFoils,1)
            if j == 1
                str = 'f';
            else
                str = 'v';
            end
            
            mazeObjects{k,j} = sprintf('%s%02d',str,nFoils(k,j));
        end
    end
    
    % Shuffle objects
    labelOrder = cell(40,1);
    % Place attended
    idx = shuffle(1:length(mazeObjects))';
    labelOrder(find(objOrder(:,4) == 1)) = mazeObjects(idx,find(attend));
    % Place unattended
    idx = shuffle(1:length(mazeObjects))';
    labelOrder(find(objOrder(:,4) == 0)) = mazeObjects(idx,find(~attend));
    
    % Write to struct
    maze(i).objects.used = mazeObjects;
    maze(i).objects.foils = foilObjects;
    maze(i).objects.order = labelOrder;
    
    
    
    %% Write Out Design
    % Unity output
    fn1 = [unityDir sprintf('obj_data_sub%g_maze%g.txt',sub,i)];
    fid1 = fopen(fn1,'w');
    
    % Text backup
    fn2 = [behDir sprintf('Configs/sub%g_config_maze%g.txt',sub,i)];
    fid2 = fopen(fn2,'w');
    fprintf(fid2,'DP junc turn attend pos x y z roty obj\n');
    
    for j = 1:numel(mazeObjects)
        % Cell array
        maze(i).design(j).DP = objOrder(j,1);
        maze(i).design(j).junc = objOrder(j,2);
        maze(i).design(j).turn = objOrder(j,3);
        maze(i).design(j).attend = objOrder(j,4);
        maze(i).design(j).pos = objOrder(j,5);
        maze(i).design(j).x = objOrder(j,6);
        maze(i).design(j).y = objOrder(j,7);
        maze(i).design(j).z = objOrder(j,8);
        maze(i).design(j).roty = objOrder(j,9);
        maze(i).design(j).label = labelOrder(j);
        
        % Text file for Unity
        fprintf(fid1,'%g,%g,%g,%g,%s\n', ...
                objOrder(j,6),objOrder(j,7),objOrder(j,8), ...
                objOrder(j,9),labelOrder{j});
    
        % Backup text file
        fprintf(fid2,'%03d %03d %+g %g %+g %+06.2f %g %+06.2f %+04.0f %s\n', ...
                objOrder(j,1),objOrder(j,2),objOrder(j,3),objOrder(j,4), ...
                objOrder(j,5),objOrder(j,6),objOrder(j,7),objOrder(j,8), ...
                objOrder(j,9),labelOrder{j});
    end
    
    labels{i} = labelOrder;
end

%% Index Objects
% Master list of object names
object = [];
for i = 1:96
    if i <= 48
        str = 'f';
        d = 0;
    else
        str = 'v';
        d = 48;
    end
    
    object{i,1} = [str sprintf('%02d',i-d)];
end

index.columns = {'Foils','Experimental','Maze','DP','Turn','Attend','Position'};
index.conditions = zeros(96,7);
dStruct = [maze.design];
for i = 1:length(object)
    
    % Object name
    index.names{i,1} = object{i};
    
    % If the object is in the maze:
    idx = find(strcmp([dStruct.label],object{i}));
    if idx
        
        % Experimental?
        index.conditions(i,2) = ~isnan(dStruct(idx).DP);
        
        % Which maze?
        if idx <= 40
            tmp = 1;
        else
            tmp = 2;
        end
        index.conditions(i,3) = tmp;
        
        % DP?
        if isnan(dStruct(idx).junc)
            tmp = 0;
        elseif dStruct(idx).junc == 0
            tmp = 1;
        else
            tmp = dStruct(idx).junc;
        end
        index.conditions(i,4) = tmp;
        
        % Turn
        index.conditions(i,5) = dStruct(idx).turn;
        
        % Attend
        index.conditions(i,6) = dStruct(idx).attend;
        
        % Position
        index.conditions(i,7) = dStruct(idx).pos;
        
    else
        
        % Its a foil
        index.conditions(i,1) = 1;
        
    end
end

    


% Save mat file
params.sub    = sub;
params.attend = attend;
params.posOff = posOff;
params.angOff = angOff;
params.DP     = DP;
params.junc   = junc;
params.turn   = turn;
params.attend = attend;
params.pos    = pos;
save([behDir sprintf('Configs/sub%02d_config.mat',sub)],'maze','params','index');






























































if 1==2
    keyboard

    % Index foils
    objects.foils = zeros(96,1);
    foilObjs = [maze(1).objects.foils; maze(2).objects.foils];
    for i = 1:length(foilObjs)
        objects.foils(find(strcmp(object,foilObjs{i})),1) = 1;
    end


    % Index exp
    objects.experiment = zeros(96,1);
    expObjs = [labels{1}(maze(1).idx.experimental==1); ...
               labels{2}(maze(1).idx.experimental==1)];
    for i = 1:length(expObjs)
        objects.experiment(find(strcmp(object,expObjs{i})),1) = 1;
    end

    % Nonexperimental
    objects.null = zeros(96,1);
    nullObjs = [labels{1}(maze(1).idx.experimental==0); ...
                labels{2}(maze(1).idx.experimental==0)];
    for i = 1:length(nullObjs)
        objects.null(find(strcmp(object,nullObjs{i})),1) = 1;
    end

    % Index maze
    objects.maze2 = zeros(96,1);
    for i = 1:length(labels{2})
        objects.maze2(find(strcmp(object,labels{2}{i})),1) = 1;
    end

    % Index turns
    turnIdx = [0 repmat([-1 -1 1 1],1,4) -1 -1 0 0 -1 -1 repmat([1 1 -1 ...
                        -1],1,4) 0]';
    turnIdx = [turnIdx; turnIdx];
    a = [labels{1}; labels{2}];
    objects.turns = zeros(96,1);
    for i = 1:length(a)
        objects.turns(find(strcmp(object,a{i})),1) = turnIdx(i);
    end

    % Index DP
    DPidx = [juncPatt(:,1); juncPatt(:,2)];
    DPidx(DPidx == 0) = 1;
    objects.DP = zeros(96,1);
    for i = 1:length(DPidx)
        objects.DP(find(strcmp(object,expObjs{i})),1) = DPidx(i);
    end

    % Index attend
    objects.attend = zeros(96,1);
    tmp = {'f','v'};
    attended = tmp(attend==1);
    for i = 1:length(expObjs)
        if strcmp(expObjs{i}(1),attended)
            objects.attended(find(strcmp(object,expObjs{i})),1) = 1;
        end
    end

    % Index visual position

    objects.pos = zeros(96,1);


    objects.labels = object;
    objects.idxMat = [objects.foils objects.experiment objects.null objects.maze2 ...
                      objects.turns];

    objects = [];
end