function fitPsycho(subs,dataType)

% Fit psychometric function to yes/no responses to dot motion stim.
clear all; close all;
if nargin < 2 error('Not enough input arguments...'); end


%% Directories
if strcmp(dataType,'afc')
    str = 'AFCFitData/';
    dType = 0;
elseif strcmp(dataType, 'yn')
    str = 'PsychoFitData/';
    dType = 1;
else
    error('Must specify ''yesno'' or ''AFC'' for data type...');
end
dataDir = ['~/Documents/Classwork/Experimental_Methods/ThresholdEst/' ...
           str];


%% Load and concatenate
catData = [];
for i = 1:length(subs)
    % Load data
    tmp = importdata([dataDir 'sub_' mat2str(subs(i)) '.dat']);
    data = [repmat(subs(i),length(tmp.data),1) tmp.data(:,2:end)];

    % Concatenate
    catData = [catData; data];
end

% Coherence values
coherences = unique(catData(:,4));
fitCoherences = linspace(min(coherences),max(coherences),100);



%% FITTING TIME
% Using the PALAMEDES toolbox cos the other one won't work...

% Set parameters for each data type:
if dType
    PF = @PAL_CumulativeNormal;
    respCol = 6;
    allResps = catData(:,respCol) > 0;
    titleStr = 'Yes/No Fits (Cum. Norm.)';
    yStr = 'p(Responding Up)';
    scaling = 0;
else
    PF = @PAL_Gumbel;
    respCol = 7;
    allResps = catData(:,respCol) == catData(:,6);
    titleStr = '2AFC Fits (Weibull)';
    yStr = 'p(Correct)';
    scaling = median(fitCoherences);
end

% Fitting parameters and options
%options = optimset('fminsearch');   % Type help optimset
%options.TolFun = 1e-09;             % Increase required precision on LL
%options.Display = 'off';            % Suppress fminsearch messages
paramsFree = [1 1 0 0];
paramsValues0 = [mean(coherences) 1/((max(coherences')-min(coherences'))/4) ...
                 0 0];
options = PAL_minimize('options');
lapseLimits = [0 1];                 % Limit range for lambda

% Fit for each sound direction:
colors = {'b','g','r'};
soundDirs = unique(catData(:,5));
hold on
for i = 1:length(soundDirs)
    % Select subset:
    tmp = [];
    tmp = catData(catData(:,5) == soundDirs(i), :);
    resps = allResps(catData(:,5) == soundDirs(i));
    realData{i} = grpstats(resps, tmp(:,4));
    
    % Group stats
    fitData{i} = [grpstats(resps, tmp(:,4), 'sum') ...
                  grpstats(resps, tmp(:,4), 'numel')];
    
    % Fitting
    [paramsValues{i}] = PAL_PFML_Fit(...
        coherences,fitData{i}(:,1),fitData{i}(:,2), ...
        paramsValues0,paramsFree,PF,'searchOptions',options, ...
        'lapseLimits',lapseLimits);
    respFit = PF(paramsValues{i}, fitCoherences');
    
    % Plot data and fit
    plot(coherences, realData{i},'.','Color',colors{i}, ...
         'MarkerSize',20);
    h(i) = plot(fitCoherences, respFit, '-', 'Color', colors{i});
    title(titleStr);
    xlabel('Coherences');
    ylabel(yStr);
    
    % Print parameters of interest
    pse(i) = PAL_CumulativeNormal(paramsValues{i}, .5, 'inverse');
    thresh(i) = PAL_CumulativeNormal(paramsValues{i}, .75, 'inverse') ...
        - PAL_CumulativeNormal(paramsValues{i}, .5, 'inverse');
    fprintf('Sound Direction %+d (Palamedes Fit):\n', soundDirs(i));
    fprintf('\tPSE: %+4.03f\n\tThresh: %4.03f\n', ...
            pse(i), thresh(i));
    
end
legend(h, 'Down', 'Flat', 'Up', 'Location', 'SE');
legend('boxoff');
hold off

return


    
    
    
    



    

