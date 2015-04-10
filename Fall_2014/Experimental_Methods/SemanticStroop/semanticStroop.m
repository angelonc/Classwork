function [resp RT correct] = semanticStroop(sub)

commandwindow


%% I/0
if nargin < 1 error('Must specify subject number'); end
mainDir = '~/Documents/Classwork/Experimental_Methods/SemanticStroop';
datDir = sprintf('%s/data/%02d', mainDir, sub);
if ~exist(datDir, 'dir') mkdir(datDir); end
str = datestr(now,'yymmdd_HHMMSS');
fn = [datDir sprintf('/sub%02d',sub) '_' str];
fid = fopen([fn '.dat'], 'w');
fprintf(fid, 'subject block trial semantic_colorIdx wordIdx display_colorIdx response correct RT\n');


%% Stimuli
% CIE Colors
colorsCIE = [.622 .346 13.8;   % red
             .538 .413 37.2;   % orange
             .423 .505 65.5;   % yellow
             .280 .617 42.3;   % green
             .152 .076 6.90;   % blue
             .387 .211 10.35]; % purple (average of red and green for now)
colorsRGB = [255 0 0; 255 127 0; 255 255 0; 0 255 0; 0 0 255; 143 0 255];
% cal = LoadCalFile('calname');
stim.colorsCIE = colorsCIE;
stim.colorsRGB = colorsRGB;

% Words
words = {'CHERRY' 'RUBY' 'BLOOD' 'RED';
         'TIGER' 'CARROTS' 'TANGERINE' 'ORANGE';
         'BANANA' 'SUN' 'LEMON' 'YELLOW';
         'LIME' 'GRASS' 'EMERALD' 'GREEN';
         'DENIM' 'OCEAN' 'SKY' 'BLUE';
         'PLUM' 'EGGPLANT' 'RAISINS' 'PURPLE'};

practice_words = {'RED' 'ORANGE' 'YELLOW' 'GREEN' 'BLUE' 'PURPLE'};
stim.words = words;
stim.practWords = practice_words;

% Design matrix
cnt = 1;
for i = 1:size(words,1)
    for j = 1:size(words,2)
        for k = 1:length(colorsCIE)
            designMat(cnt,1) = i;
            designMat(cnt,2) = j;
            designMat(cnt,3) = k;
            cnt = cnt + 1;
        end
    end
end

% Practice Mat
cnt = 1;
for i = 1:length(practice_words)
    for j = 1:length(colorsRGB)
        practMat(cnt,1) = i;
        practMat(cnt,2) = j;
        cnt = cnt + 1;
    end
end

% Parameters
params.nBlocks = 3;
params.nTrials = length(designMat);
params.nPractTrials = length(practMat);
params.respWait = 3;
params.ITI = .5;
params.keysOfInterest = zeros(1,256);
params.keysOfInterest(KbName({'1!','2@','3#','4$','5%','6^', 'ESCAPE'})) = 1;

save([fn '_params.mat'], 'params', 'stim');


%% Start Experiment
try
    % Open window
    screen = max(Screen('Screens'));
    black = BlackIndex(screen);
    [win, rect] = Screen('OpenWindow', screen, [], [0 0 600 400]);
    Screen('TextFont', win, 'Arial');
    Screen('TextSize', win, 24);
    HideCursor;
    
    
    % Display instructions
    respStr = '"1" = RED\n"2" = ORANGE\n"3" = YELLOW\n"4" = GREEN\n"5" = BLUE\n"6" = PURPLE\n\n';
    instructions1 = [sprintf('INSTRUCTIONS (Part 1):\n\nOn each trial you will see different words of different colors.\nYour task is simply to report the color on the screen of the word presented.\n\nYou will report word color by pressing keys on the keyboard, as follows:\n\n') sprintf(respStr) '(press any key to continue to part 2)'];
    Screen('FillRect', win, black, rect);
    DrawFormattedText(win, instructions1, 'center', 'center', [255 255 255]);
    Screen('Flip', win);
    WaitSecs(.5);
    KbWait;
    
    instructions2 = sprintf('INSTRUCTIONS (Part 2):\n\nBefore you begin the real task, lets practice making correct color responses.\n\n This practice will show words that correspond to each color.\nYour job is to press the key that corresponds to the color of the word presented\n(eg. of the color of the word is red, press 1; if it is orange, press 2, and so on)\n\nPress any key to continue to the practice.');
    Screen('FillRect', win, black, rect);
    DrawFormattedText(win, instructions2, 'center', 'center', [255 255 255]);
    Screen('Flip', win);
    WaitSecs(.5);
    [~,keyCode] = KbWait;
    skipPract = any(find(keyCode == 1) == KbName('ESCAPE'));
    WaitSecs(1);
    
    
    %% Practice
    if ~skipPract
        % Restrict keys/start queue
        KbQueueCreate(-1, params.keysOfInterest);
        
        % Get ready
        practTrials = params.nPractTrials;
        idx = Shuffle(1:practTrials);
        resp = zeros(1,practTrials);
        RT = zeros(1,practTrials);
        corr = zeros(1,practTrials);
        startTime = GetSecs;
        
        % GO!
        for i = 1:practTrials
            % Write conditions to file
            fprintf(fid, '%02d %02d %03d %g %g %g ', ...
                    sub, 0, i, practMat(i,1), 0, practMat(i,2));
            
            KbQueueStart;
            Screen('FillRect', win, black, rect);
            DrawFormattedText(win, practice_words{practMat(idx(i),1)}, ...
                              'center', 'center', ...
                              colorsRGB(practMat(idx(i),2),:));
            ts = Screen('Flip', win);
            
            % Go to next trial after specified wait time
            dt = 0;
            pressed = 0;
            while 1
                dt = GetSecs - ts;
                
                % Check for a press
                [pressed, pressT] = KbQueueCheck;
                
                if dt >= params.respWait
                    % If time is up, break out of the loop
                    resp(i) = 1;
                    RT(i) = NaN;
                    break;
                elseif pressed                   
                    % If a key is pressed:
                    disp('Key pressed');
                    if any(find(pressT > 0) == KbName('ESCAPE'));
                        DrawFormattedText(win, 'Quitting...', 'center', 'center', [255 255 255]);
                        Screen('Flip',win);
                        WaitSecs(1);
                        Screen('CloseAll');
                        return;
                    end
                    resp(i) = find(pressT > 0);
                    RT(i) = pressT(pressT > 0) - ts;
                    break;
                end
                WaitSecs(.001);
            end
            
            % Flush the queue
            KbQueueStop;
            KbQueueFlush;
            
            % Show feedback
            respKey = KbName(resp(i));
            respKey = str2double(respKey(1));
            correct(i) = respKey == practMat(idx(i),2);
            if correct(i)
                feedColor = [0 255 0];
                feedStr = 'CORRECT';
            else
                feedColor = [255 0 0];
                feedStr = 'INCORRECT';
            end
            DrawFormattedText(win, feedStr, 'center', 'center', feedColor);
            Screen('Flip', win);
            WaitSecs(.5);
            
            % Write response to file
            fprintf(fid, '%g %g %.4f\n', respKey, correct(i), RT(i));
            
            % Fill with black till next trial starts
            Screen('FillRect', win, black, rect);
            Screen('Flip', win);
            WaitSecs(params.ITI);
        end
        KbQueueRelease;
        
        % Display percent correct and feedback
        practEndStr = [sprintf('Nice job! You got %02.2f%% of the trials correct!\n\n', mean(correct)*100) ...
                       sprintf('The real task will be next. ')];
    else
        practEndStr = [];
    end
    
    % Put up instructions for real task
    instructions3 = sprintf('In this task there will be %g blocks with %03d trials in each.\n\nPress any key when you are ready to begin block 1.', params.nBlocks, params.nTrials);
    DrawFormattedText(win, [practEndStr instructions3], 'center', 'center', [255 255 255]);
    Screen('Flip', win);
    WaitSecs(.5);
    KbWait;
    Screen('FillRect', win, black, rect);
    Screen('Flip', win);
    WaitSecs(1);
    
    %% Real Task
    % For each block
    for b = 1:params.nBlocks
        % Pause before blocks
        if b > 1
            str = sprintf('You are about to begin block %g/%g. Press any key to start.', b, params.nBlocks);
            DrawFormattedText(win, str, 'center', 'center', [255 255 255]);
            Screen('Flip', win);
            WaitSecs(.5);
            KbWait;
            Screen('FillRect', win, black, rect);
            Screen('Flip', win);
            WaitSecs(1);
        end
        
        % Start the queue
        KbQueueCreate(-1, params.keysOfInterest);
        
        % Shuffle the design matrix
        idx = shuffle(designMat);
        
        % For each trial in a block
        for t = 1:params.nTrials
            % Write conditions to file
            fprintf(fid, '%02d %02d %03d %g %g %g ', ...
                    sub, b, t, idx(t,1), idx(t,2), idx(t,3));
            
            KbQueueStart;
            Screen('FillRect', win, black, rect);
            DrawFormattedText(win, words{idx(t,1),idx(t,2)}, ...
                              'center', 'center', ...
                              colorsRGB(idx(t,3),:));
            ts = Screen('Flip', win);
            
            % Go to next trial after specified wait time
            dt = 0;
            pressed = 0;
            while 1
                dt = GetSecs - ts;
                
                % Check for a press
                [pressed, pressT] = KbQueueCheck;
                
                if dt >= params.respWait
                    % If time is up, break out of the loop
                    resp(t) = 1;
                    RT(t) = NaN;
                    break;
                elseif pressed                   
                    % If a key is pressed:
                    disp('Key pressed');
                    if any(find(pressT > 0) == KbName('ESCAPE'));
                        fprintf(fid, '\nUser quit.\n');
                        fclose(fid);
                        DrawFormattedText(win, 'Quitting...', 'center', 'center', [255 255 255]);
                        Screen('Flip',win);
                        WaitSecs(1);
                        Screen('CloseAll');
                        return;
                    end
                    resp(t) = find(pressT > 0);
                    RT(t) = pressT(pressT > 0) - ts;
                    break;
                end
                WaitSecs(.001);
            end
            
            % Flush the queue
            KbQueueStop;
            KbQueueFlush;
            
            % Code responses and calculate correct
            respKey = KbName(resp(t));
            respKey = str2double(respKey(1));
            correct(t) = respKey == idx(t,3);
            
            % Write to file
            fprintf(fid, '%g %g %.4f\n', respKey, correct(t), RT(t));                        
            
            % Fill with black till next trial starts
            Screen('FillRect', win, black, rect);
            Screen('Flip', win);
            WaitSecs(params.ITI);
        end
        KbQueueRelease;
    end
    
    
catch err
    rethrow(err)
    sca
    fclose('all');
    KbQueueRelease;
    keyboard;
end
fclose(fid);
sca

return;
sca











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out, idx] = shuffle(in)

% Shuffle along longest dimension
len = length(in);
idx = randperm(len)';
out = in(idx,:);

return










