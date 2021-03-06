function [location, re_loc, post, re_loc_ncc post_nocc, pv_loc] = place_recon(mouse,CC)

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
    lower = 0;
    upper = 35;
    minV = 3;
    bins = 0:.5:35;
    nBins = length(bins)-1;
    kernelSize = 3;
    mask = makeMask(0);
else
    % Human data
    load delivery_subject1_session1.mat
    
    % Extract relevant info
    x = [sessionData.behav.x];
    y = [sessionData.behav.y];
    sr = 1/sessionData.epochSize;
    
    % Calculate velocity
    v = abs([0 diff(sqrt(x.^2 + y.^2))*1000]); % vu/s
    
    % Get spikes
    spikes = sessionData.fr';
    
    % Make bins
    upper = ceil(max([max(x(:)) max(y(:))]));
    lower = floor(min([min(x(:)) min(y(:))]));
    nBins = 70;
    bins = linspace(lower,upper,nBins + 1);
    
    % Params
    minV = .001;
    kernelSize = 12;
    mask = ones(nBins,nBins);
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


%% MAPS
% Index the training period
train = .75;

train_idx = zeros(1,length(x));
train_idx(1:round(length(x)*train)) = 1;
train_idx = logical(train_idx);

v_train = v(train_idx & move_idx);
x_train = x(train_idx & move_idx);
y_train = y(train_idx & move_idx);
spikes_train = spikes(:,train_idx & move_idx);

% Make sampling maps
[N,C] = hist3([x_train' y_train'],'Edges',edges);
s_map = N(1:end-1,1:end-1) / sr; % Convert to seconds
p_map = s_map / sum(s_map(:)) .* mask;
s_map = imfilter(s_map,gKern,'same');
s_map = s_map .* mask;

% Activity maps
disp('Activity maps...');

nCells = size(spikes,1);
spike_map = zeros(nBins,nBins,nCells);
f_map = zeros(nBins,nBins,nCells);
pf_map = zeros(nBins,nBins,nCells);
for i = 1:size(spikes,1)
    % Make an activity map
    spike_map(:,:,i) = genFiringFields(x_train,y_train,...
                                       spikes_train(i,:),bins);
    f_map(:,:,i) = imfilter(spike_map(:,:,i),gKern);
    pf_map(:,:,i) = f_map(:,:,i) ./ s_map;
end


%% RECONSTRUCTIONS
disp('Reconstruction...');

timeBin = 1;
samps_per_bin = round(sr * timeBin);
n_bins = floor(length(x(~train_idx))/4);

x_test = x(~train_idx); 
y_test = y(~train_idx);
v_test = v(~train_idx);
spikes_test = spikes(:,~train_idx);


n = zeros(1,1,size(spikes,1));
t_bin = zeros(n_bins,samps_per_bin);
live_plot = 1;
pv_loc = [];
% For each time bin
for t = 1:n_bins
    disp(['Time ' num2str(t)]);
    t_bin(t,:) = [((t-1) * samps_per_bin) + 1 :...
            ((t-1) * samps_per_bin) + samps_per_bin];
    
    % Current location
    loc(t,:) = [mean(x_test(t_bin(t,:))) mean(y_test(t_bin(t,:)))];
    
    
    %% BAYESIAN
    % Continuity constraint around previous position
    if t == 1
        CC = ones(70,70);
    else
        % Centered on last position
        %Mu = fliplr([mean(x_test(t_bin(t-1,:))) mean(y_test(t_bin(t-1,:)))]);
        Mu = fliplr(re_loc(t-1,:));
        pos = zeros(70,70);
        pos(Mu(2),Mu(1)) = 1;
        % Width is present velocity
        Sigma = abs(mean(v_test(t_bin(t,:))))/params.binSize;
        
        kern = fspecial('gaussian',7*[ceil(Sigma) ceil(Sigma)],Sigma);
        CC = imfilter(pos,kern);

%          Sigma = repmat(Sigma,1,2);
%          CC = reshape(mvnpdf([X1(:) Y1(:)],Mu,Sigma),71,71);
%          CC = CC(1:end-1,1:end-1) ./ sum(CC(:));
    end
   
    for unit = 1:size(spikes,1)
        % Make population vector
        n(unit) = sum(spikes_test(unit,t_bin(t,:)));
        pf_power(:,:,unit) = (pf_map(:,:,unit) .^ n(unit));
    end
    
    post_nocc{t} = p_map .* prod(pf_power,3) .* exp(-params.timeBin .* sum(pf_map,3));
    [~,argmax_idx] = max(post_nocc{t}(:));
    [re_loc_ncc(t,1),re_loc_ncc(t,2)] = ind2sub(size(post_nocc{t}),argmax_idx);
    
    % Calculate posterior
    post{t} = CC .* p_map .* prod(pf_power,3) .* exp(-params.timeBin .* sum(pf_map,3));
    [~,argmax_idx] = max(post{t}(:));
    [re_loc(t,1),re_loc(t,2)] = ind2sub(size(post{t}),argmax_idx);
        
    
    %% POP VECTOR
    warning off

    for ii = 1:length(bins)-1
        for jj = 1:length(bins) - 1
            dist_map(ii,jj) = pdist([n(:)'; squeeze(spike_map(ii,jj,:))'],'cosine');
        end
    end
    [~,argmax_idx] = min(dist_map(:));
    [pv_loc(t,1),pv_loc(t,2)] = ind2sub(size(dist_map),argmax_idx);

    
    if 0
        ax = subplot(1,4,1);
        plot(loc(t,2),loc(t,1),'o');
        xlim([lower upper])
        ylim([lower upper])
        set(ax,'YDir','reverse');
        subplot(1,4,2)
        imagesc(CC)
        subplot(1,4,3)
        imagesc(post{t})
        subplot(1,4,4)
        imagesc(1 - dist_map)
        keyboard
    end
    
end



























































%% SUBFUNCS


function [f_map] = genFiringFields(x,y,spikes,bins)

n_bins = length(bins);
f_map = zeros(n_bins-1,n_bins-1);
for xx = 1:n_bins - 1
    for yy = 1:n_bins - 1
        x_bin = [bins(xx) bins(xx+1)];
        y_bin = [bins(yy) bins(yy+1)];
        
        t_idx = (x > x_bin(1) & x <= x_bin(2)) & ...
                (y > y_bin(1) & y <= y_bin(2));
        
        f_map(xx,yy) = sum(spikes(t_idx));
    end
end



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




function [mask] = makeMask(vals)

cx = 35.5; cy = 35.5;
r = 35.5;
ix = 71;
iy = 71;
[mx,my]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
c_mask=((mx.^2+my.^2)<=r^2);
mask = +c_mask(1:end-1,1:end-1);
mask(mask == 0) = vals;



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    