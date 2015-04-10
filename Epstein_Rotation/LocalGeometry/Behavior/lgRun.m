function lgRun(sub,block,trial)

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
% Directories
Dir = ['~/Documents/Classwork/Epstein_Rotation/LocalGeometry/' ...
       'Behavior/'];
datDir = [Dir sprintf('Data/%02d/',sub)];

% Load in images and order
fprintf('Loading images...\n');
load([Dir 'Configs/cropped_ims.mat']);
fprintf('Loading order...\n');
load([datDir 'test_order.mat']);




%% SETUP
% Open a new screen
screen = max(Screen('Screens'));
[win, winRect] = Screen('OpenWindow',screen,[bg bg bg],[0 0 800 600]);
IFI = Screen('GetFlipInterval', win);
center = winRect(3:4)./2;
Screen('BlendFunction',win);

% Make the textures
makeTextures(win,winRect,im,im_files,rect,center,bg)




%% EXPERIMENT
% Parameters
params.nBlocks   = length(beh.order);
params.nTrials   = length(beh.order{1});
params.imTime    = .5;
params.fixTime   = 12;
params.resp1Time = NaN;
params.resp2Time = NaN;
params.ITI       = 1;
params.keysOfInterest(KbName({'1!','2@','ESCAPE'})) = 1;

% For each block
for i = block:params.nBlocks
    % Display instructions
    if i == 1
        DrawFormattedText(win,'Instructions here','center','center',[0 ...
                            0 0]);
        Screen('Flip',win);
    else
        DrawFormattedText(win,'Break','center','center',[0 ...
                            0 0]);
        Screen('Flip',win);
    end
    
    % Setup KbQueue
    KbQueueCreate(-1,params.keysOfInterest);
    
    % Wait for trigger (or keypress)
    if ~params.fMRI
        ts.start = KbWait;
    else
        % Get trigger
    end
    
    % For each trial
    for j = trial:params.nTrials
        % Start the queue
        KbQueueStart;
        
        % Draw texture
        idx = beh.order{i}(j,:);
        Screen('DrawTexture',win, ...
               texture(idx(1),idx(2),idx(3)),[], ...
               tex_rect{idx(1),idx(2),idx(3)});
        drawFixation(win,bg,center);
        ts.flip(i) = Screen('Flip',win);
        
        % Do you recognize it?
        DrawFormattedText(win,'Did you see this object in the maze?\n[1] = YES\n[2] = NO ');
        ts.q1(i) = Screen('Flip',win, ts.flip(i) + params.imTime);
        
        dt = 0;
        pressed = 0;
        while 1
            dt = GetSecs - ts;
            
            [pressed, ts.r1] = KbQueueCheck;
            
            
            WaitSecs(.0001);
        end
    end
end


Screen('CloseAll');
Screen('Close');
keyboard








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function drawFixation(win,bg,center)

Screen('DrawDots',win,center,10,[255 255 255],[0 0],2);
Screen('DrawDots',win,center,7,[0 0 0],[0 0],2);
Screen('Flip',win);



function clearScreen(win,bg,winRect)

Screen('FillRect',win,[bg bg bg],winRect);
Screen('Flip',win);



function makeTextures(win,winRect,im,im_files,rect,center,bg)

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
clearScreen(win,bg,winRect);
