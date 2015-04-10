function run_localgeom(sub,block,trial)

dbstop if error

%% Variables
if ~exist('sub','var')
    sub = 99;
end
if ~exist('block','var')
    block = 1;
end
if ~exist('trial','var')
    trial = 1;
end

commandwindow




%% IO
try
    % Directories
    Dir = ['~/Documents/Classwork/Epstein_Rotation/LocalGeometry/' ...
           'Behavior/'];
    datDir = [Dir sprintf('Data/%02d/',sub)];
    if ~exist(datDir,'dir') error('ERROR: Make sure this participant has a design file made!'); end;

    % Load in images and order
    fprintf('Loading images...\n');
    load([Dir 'Configs/cropped_ims.mat']);
    fprintf('Loading order...\n');
    load([datDir 'test_order.mat']);
    fprintf('Loading design...\n');
    load([Dir sprintf('Configs/sub%02d_config.mat',sub)]);
    
    % Open a file
    today = datestr(now,'yymmdd_HHMMSS');
    fn = [datDir sprintf('sub_%02d_%s',sub,today)];
    fid = fopen([fn '.dat'], 'w');
    fprintf(fid,['sub block trial obj foil exp maze DP ' ...
                 'turn attend pos resp1 resp2 resp3 RT1 RT2 RT3 correct1 correct2 correct3\n']);



    %% SETUP
    % Open a new screen
    fprintf('Starting PTB...\n');
    screen = max(Screen('Screens'));
    [win, winRect] = Screen('OpenWindow',screen,[bg bg bg],[0 0 1000 800]);
    IFI = Screen('GetFlipInterval', win);
    center = winRect(3:4)./2;
    Screen('BlendFunction',win);
    Screen(win,'TextSize',20);

    % Make the textures
    [texture, tex_rect] = makeTextures(win,winRect,ims,im_files,rect,center,bg);
    clearScreen(win,bg,winRect);
    Screen('Flip',win);
    
    % Instruction text
    instructions{1} = sprintf(['You will be presented with images of objects you have' ...
              ' encountered\nor have not encountered in the museums.\n\n']);
    instructions{2} = sprintf(['For each object, you will indicate:\n'...
              '\t1) Whether you have seen it in a museum [1] or not [2].\n'...
              'If you have seen it:\n'...
              '\ta) Did you see it in the brick [1] or the wood museum [2]?\n'...
              '\tb) Did you make a left [1] or right [2] turn at that object?\n\n']);
    instructions{3} = sprintf(['Please try to answer as quickly and accurately as possible!\n\n'...
              'Press any key to start.']);




    %% EXPERIMENT

    % Parameters
    params.nBlocks     = length(beh.order);
    params.nTrials     = length(beh.order{1});
    params.imTime      = 1;
    params.fixTime     = 12;
    params.resp1Time   = NaN;
    params.resp2Time   = NaN;
    params.ITI         = .5;
    params.queryHeight = center(2) + 310;

    % Restrict keys
    params.keysOfInterest = zeros(1,256);
    params.keysOfInterest(KbName({'1!','2@'})) = 1;



    %% For each block
    for j = block:params.nBlocks
        clear resp correct
        
        %% Block Setup
        if j == 1
            drawInstructions(win, instructions, center);
            Screen('Flip',win);
        else
            str = sprintf('Feel free to take a break.\n\nPress any key to start block %g/%g',...
                          j,params.nBlocks);
            DrawFormattedText(win,str,'center','center',[0 ...
                                0 0]);
            Screen('Flip',win);
        end   
        WaitSecs(1);
        
        % Wait for trigger (or keypress)
        ts.start = KbWait;
        
        % Setup KbQueue
        KbQueueCreate([],params.keysOfInterest);
        
        
        
        %% For each trial
        for i = trial:params.nTrials
            
            % Get trial start time
            ts.trialStart(i) = GetSecs;
 
            %% Draw texture
            idx = beh.order{j}(i,:);
            Screen('DrawTexture',win, ...
                   texture(idx(1),idx(2),idx(3)),[], ...
                   tex_rect{idx(1),idx(2),idx(3)});
            drawFixation(win,bg,center);
            DrawFormattedText(win,['Have you seen this object?\n[1] ' ...
                                '= YES\n[2] = NO'], 'center',params.queryHeight);
            ts.flip(i) = Screen('Flip',win, ts.trialStart(i) + params.ITI);

            

            %% Get responses
            % Check for responses
            keys = getResponse;


            % Record responses
            ts.query(i,1) = ts.flip(i);
            resp(i,1) = find(keys > 0) - 29;
            RT(i,1) = keys(keys > 0) - ts.query(i,1);
            
            % If you did see it, ask which maze and turn you made:
            if find(keys > 0) == KbName('1!')
                % Do you recognize it?
                Screen('DrawTexture',win, ...
                       texture(idx(1),idx(2),idx(3)),[], ...
                       tex_rect{idx(1),idx(2),idx(3)});
                DrawFormattedText(win,'Did you see it the brick or wood museum?\n[1] = BRICK MAZE\n[2] = WOOD MAZE', ...
                                  'center',params.queryHeight);
                ts.query(i,2) = Screen('Flip',win);
                
                % Get maze response
                keys = getResponse;
                
                resp(i,2) = find(keys > 0) - 29;
                RT(i,2) = keys(keys > 0) - ts.query(i,2);
                
                % What turn was made?
                Screen('DrawTexture',win, ...
                       texture(idx(1),idx(2),idx(3)),[], ...
                       tex_rect{idx(1),idx(2),idx(3)});
                DrawFormattedText(win,'What turn was made at this object?\n[1] = LEFT \n[2] = RIGHT', ...
                                  'center',params.queryHeight);
                ts.query(i,3) = Screen('Flip',win);
                
                % Get maze response
                keys = getResponse;
                
                resp(i,3) = find(keys > 0) - 29;
                RT(i,3) = keys(keys > 0) - ts.query(i,3);
            else
                ts.query(i,2) = NaN;
                ts.query(i,3) = NaN;
                
                resp(i,2) = NaN;
                resp(i,3) = NaN;
                
                RT(i,2) = 0;
                RT(i,3) = 0;
            end
            
            % Current index
            objIdx = find(strcmp(index.names,beh.names{i,j}));
            
            % Write to file
            fprintf(fid,'%02d %g %02d %s %g %g %g %g %+g %g %+g %g %g %g %07.4f %07.4f %07.4f ',...
                    sub,...
                    i,...
                    j,...
                    beh.names{i,j},...
                    index.conditions(objIdx,1),...
                    index.conditions(objIdx,2),...
                    index.conditions(objIdx,3),...
                    index.conditions(objIdx,4),...
                    index.conditions(objIdx,5),...
                    index.conditions(objIdx,6),...
                    index.conditions(objIdx,7),...
                    resp(i,1),...
                    resp(i,2),...
                    resp(i,3),...
                    RT(i,1),...
                    RT(i,2),...
                    RT(i,3));
                    
           
            
            %% Correct responses?
            
            
            % Recognition?
            mazeObjs = index.conditions(:,4) > 0;
            if resp(i,1) == 1 & mazeObjs(objIdx) == 1
                correct(i,1) = 1;
            elseif resp(i,1) == 2 & mazeObjs(objIdx) == 0
                correct(i,1) = 1;
            else
                correct(i,1) = 0;
            end
            
            % Maze?
            if ~isnan(resp(i,2))
                correct(i,2) = index.conditions(objIdx,3) == resp(i,2);
            else
                correct(i,2) = NaN;
            end
            
            % Turn?
            if resp(i,3) == 1 & index.conditions(objIdx,5) == -1
                correct(i,3) = 1;
            elseif resp(i,3) == 2 & index.conditions(objIdx,5) == 1
                correct(i,3) = 1;
            elseif isnan(resp(i,3))
                correct(i,3) = NaN;
            else
                correct(i,3) = 0;
            end
            
            % Write correct to file
            fprintf(fid,'%g %g %g\n', correct(i,1), correct(i,2), correct(i,3));
            
            % Draw Fixation
            clearScreen(win,bg,winRect);
            drawFixation(win,bg,center);
            Screen('Flip',win);
        end
        KbQueueRelease();
    end

catch err
    rethrow(err);
    KbQueueRelease();
    Screen('CloseAll');
    Screen('Close');
    keyboard
end

Screen('CloseAll');
Screen('Close');
ShowCursor;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function drawFixation(win,bg,center)

Screen('DrawDots',win,center,10,[255 255 255],[0 0],2);
Screen('DrawDots',win,center,7,[0 0 0],[0 0],2);



function ts = clearScreen(win,bg,winRect)

Screen('FillRect',win,[bg bg bg],winRect);



function [texture, tex_rect] = makeTextures(win,winRect,ims,im_files,rect,center,bg)

% Make textures
DrawFormattedText(win,'Loading textures...','center','center', [0 0 ...
                    0]);
Screen('Flip',win);
fprintf('Making textures...\n');
texture = nan(2,48,7);
for i = 1:length(ims)
    % First filenames by id number, type, angle
    str = strread(im_files(i).name,'%s','delimiter','_*.');
    type = str{1}(1);
    num = str2num(str{1}(2:3));
    angle = str2num(str{2});
    if strcmp(type,'f')
        idx = 1;
    else strcmp(type,'v')
        idx = 2;
    end
    
    % Assign textures based on (number,type,angle)
    texture(idx,num,(angle/30)+1) = Screen('MakeTexture',win,ims{i});
    
    tex_rect{idx,num,(angle/30)+1} = CenterRectOnPoint(rect(i,:),center(1),center(2));
end


function keys = getResponse

KbQueueStart();

pressed = 0;
while ~pressed
    [pressed, keys] = KbQueueCheck;
    WaitSecs(.0001);
end
KbQueueStop();

function drawInstructions(win, instructions, center)

DrawFormattedText(win,instructions{1},'center',center(2)-150,[0 0 0]);
DrawFormattedText(win,instructions{2},center(1)-300,'center',[0 0 0]);
DrawFormattedText(win,instructions{3},'center',center(2)+60,[0 0 0]);




%% Unused Code

if 1==2
    flag = 0;
    while ((GetSecs - ts.flip(i)) < (params.imTime - .16))
        % Check for a press
        [pressed, keys] = KbQueueCheck;
        WaitSecs(.0001);
        
        if pressed
            break;
        end
    end
    
    % Clear stimulus
    clearScreen(win,bg,winRect);
    drawFixation(win,bg,center);
    Screen('Flip',win,ts.flip(i) + params.imTime);
    
    % Wait for response if there hasn't been one yet
    if ~pressed
        % Do you recognize it?
        DrawFormattedText(win,'Did you see this object in a maze?\n[1] =YES\n[2] = NO ', ...
                          'center','center');
        Screen('Flip',win);
        keys = getResponse;
    end
end
