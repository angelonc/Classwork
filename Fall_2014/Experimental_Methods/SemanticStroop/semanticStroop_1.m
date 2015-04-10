function [resp RT] = semanticStroop

clear all


%% I/0
% if nargin < 1 error('Must specify subject number'); end
% mainDir = '~/Documents/Classwork/Experimental_Methods/SemanticStroop/';
% datDir = [mainDir 'data/'];
% str = datestr(now,'yymmdd_HHMMSS');
% FN = [mainDir datDir sprintf('sub%02d',sub) '_' str '.dat'];


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

% Words
words = {'CHERRY' 'RUBY' 'BLOOD' 'RED';
         'TIGER' 'CARROTS' 'TANGERINE' 'ORANGE';
         'BANANA' 'SUN' 'LEMON' 'YELLOW';
         'LIME' 'GRASS' 'EMERALD' 'GREEN';
         'DENIM' 'OCEAN' 'SKY' 'BLUE';
         'PLUM' 'EGGPLANT' 'RAISINS' 'PURPLE'};

practice_words = {'SHOE';'BOAT';'HORN';'FACE';'LONG';
                  'TORN';'LEGS';'LATE';'YEAR';'BABY';
                  'HAMS';'THOR';'FLAY';'SAVE';'TURN';
                  'HARM';'WORD';'WITH'};

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

% Parameters
nBlocks = 3;
nTrials = length(designMat);
practTrials = length(practice_words);
waitTime = 3;
ITI = .5;


%% Start Experiment
try
% Open window
    screen = max(Screen('Screens'));
    black = BlackIndex(screen);
    [win, rect] = Screen('OpenWindow', screen, [], [0 0 600 400]);
    %Screen('Preference', 'DefaultFontSize');
    Screen('TextFont', win, 'Arial');
    Screen('TextSize', win, 12);
    %HideCursor;


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
keyboard;
skipPract = any(find(keyCode == 1) == KbName('ESCAPE'));
    

    %% Practice
 
        % Restrict keys/start queue
        keysOfInterest = zeros(1,256);
        keysOfInterest(KbName({'1!','2@','3#','4$','5%','6^'})) = 1;
        KbQueueCreate(-1, keysOfInterest);
        
        % Get ready
        wordIdx = shuffle(1:practTrials);
        colorIdx = shuffle(repmat(1:6,1,practTrials/6));
        resp = zeros(1,practTrials);
        RT = zeros(1,practTrials);
        corr = zeros(1,practTrials);
        startTime = GetSecs;
        
        % GO!
        for i = 1:practTrials
            KbQueueStart;
            Screen('FillRect', win, black, rect);            
            [~,~,bounds] = DrawFormattedText(win, practice_words{wordIdx(i)}, ...
                                             'center', 'center', ...
                                             colorsRGB(colorIdx(i),:));
            ts = Screen('Flip', win);
            
            % Go to next trial after specified wait time
            dt = 0;
            pressed = 0;
            while 1
                dt = GetSecs - ts;
                
                % Check for a press
                [pressed, pressT] = KbQueueCheck;
                if dt >= waitTime
                    resp(i) = 1;
                    RT(i) = NaN;
                    break;
                elseif pressed 
                    disp('Key pressed');
                    keyboard;
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
            corr(i) = str2double(respKey(1)) == colorIdx(i)
            if corr(i)
                feedColor = [0 255 0];
                feedStr = 'CORRECT';
            else
                feedColor = [255 0 0];
                feedStr = 'INCORRECT';
            end
            %DrawFormattedText(
            Screen('FrameRect', win, feedColor, bounds, 1);
            Screen('Flip', win);
            WaitSecs(.5);
            
            % Fill with black till next trial starts
            Screen('FillRect', win, black, rect);            
            Screen('Flip', win);
            WaitSecs(ITI);
        end
        KbQueueRelease;

    
% Display percent correct and feedback

catch
    sca
    KbQueueRelease;
    keyboard;
end
sca

return;








