function [location, recon, bins] = place_recon(mouse,CC)

%% BAYESIAN SPATIAL RECONSTRUCTION
% This function calculates the estimates spatial position for a
% given time interval using ensemble firing activity. We output
% several model predictions, a population vector model, a Bayes
% reconstruction (with or without continuity constraint) and a
% Bayesian update model (prior = previous prior).
%
%
%
%% function [location, recon, bins] = place_recon(mouse,CC,update)
% OUTPUTS:
%  location = vector of average location per time bin
%  recon    = struct with reconstruction vectors for each model
%  bins     = bin dimensions (for plotting)
% INPUTS:
%  mouse    = boolean for running mouse vs human recon
%  CC       = boolean for continuity constraint

if nargin < 3
    mouse = 1;
    CC = 0;
end


%% SETUP

if mouse
    % Mouse data
    load Neutral_s1_d1.mat
    
    x = pos.rawX;
    y = pos.rawY;
    ts = pos.ts;
    sr = 1000/mean(diff(ts));

    % Calculate velocity
    v = abs([0 diff(sqrt(x.^2 + y.^2))*1000]); % cm/s
    
    % Align spikes
    spikes = alignSpikeTimes(ts,unit);
    
    % Params
    minV = 3;
    bins = 0:.5:35;
    kernelSize = 3;
else
    % Human data
    load delivery_subject1_session1.mat
    
    % Extract relevant info
    x = [sessionData.behav.x];
    y = [sessionData.behav.y];
    sr = 1/sessionData.epochSize;
    
    % Calculate velocity
    v = abs([0 diff(sqrt(x.^2 + y.^2))*1000]); % vu/s
    
    spikes = sessionData.fr';
    
    % Params
    minV = .001;
    kernelSize = 12;
    
    % Make bins
    upper = ceil(max([max(x(:)) max(y(:))]));
    lower = floor(min([min(x(:)) min(y(:))]));
    nBins = 70;
    bins = linspace(lower,upper,nBins + 1);
end

% Index periods of motion
move_idx = v > minV;

% Make a meshgrid for later
[X1, Y1] = meshgrid(bins,bins);
edges{1} = bins;
edges{2} = bins;

binSize = mean(diff(bins));

% Gaussian kernel
sigma = (kernelSize / 2*sqrt(2*log(2)))/binSize;
gKern = fspecial('gaussian', 7*[ceil(sigma) ceil(sigma)], sigma);








































































function [spike_idx] = alignSpikeTimes(ts,unit)

fr_limit = .12; % Exclude cells that fire below 12Hz
                % (50spikes/10min)
dur = (ts(end) - ts(1)) / 1000; % Duration in seconds

cell_cnt = 0;
% For each cell
for cell = 1:numel(unit)
    % Only look at cells firing over .12Hz
    if length(unit{cell})/dur > fr_limit
        cell_cnt = cell_cnt + 1;
        
        % For each spike
        spike_cnt = 0;
        spike_idx(cell_cnt,:) = zeros(1,length(ts));
        for s = 1:length(unit{cell})
            [~,i] = min(abs(ts - unit{cell}(s)));
            spike_idx(cell_cnt,i) = spike_idx(cell_cnt,i) + 1;
        end
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    