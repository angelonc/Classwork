function [loc, recon, bins] = place_recon_update(mouse,live_plot)

close all
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

if nargin < 2
    mouse = 1;
    live_plot = 1;
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
    mask = makeMask(NaN);
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

if ~exist('activity_maps.mat','file')
    
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
        pf_map(isnan(pf_map)) = 0;
    end

    save('activity_maps.mat','spike_map','f_map','pf_map');

else
    load activity_maps.mat
end


%% RECONSTRUCTIONS
disp('Reconstruction...');

timeBin = 1;
samps_per_bin = round(sr * timeBin);
n_bins = floor(length(x(~train_idx))/sr);

x_test = x(~train_idx); 
y_test = y(~train_idx);
v_test = v(~train_idx);
spikes_test = spikes(:,~train_idx);


n = zeros(1,1,size(spikes,1));
t_bin = zeros(n_bins,samps_per_bin);
pv_loc = [];
% For each time bin
for t = 1:n_bins
    disp(['Time ' num2str(t)]);
    t_bin(t,:) = [((t-1) * samps_per_bin) + 1 :...
            ((t-1) * samps_per_bin) + samps_per_bin];
    
    % Current location
    loc(t,:) = [mean(x_test(t_bin(t,:))) mean(y_test(t_bin(t,:)))];
    
    
    %% BAYESIAN
    % Get current activity
    for unit = 1:size(spikes,1)
        % Make population vector
        n(unit) = sum(spikes_test(unit,t_bin(t,:)));
        pf_power(:,:,unit) = (pf_map(:,:,unit) .^ n(unit));
    end
    
    % Set prior
    if t == 1
        prior = p_map;
        prior(isnan(prior)) = 0;
        kern = gKern;
        CC = ones(70,70);
    else
        % Update prior
        prior = recon.update.posterior(:,:,t-1);
        
        % Update CC by centering on previous prediction, width
        % scales w. velocity
        Mu = recon.cc.loc(t-1,:);
        pos = zeros(70,70);
        pos(Mu(1),Mu(2)) = 1;
        Sigma = abs(mean(v_test(t_bin(t,:))))/params.binSize;
        kernCC = fspecial('gaussian',7*[ceil(Sigma) ceil(Sigma)],Sigma);
        CC = imfilter(pos,kernCC) ./ trapz(CC(:));
        
    end
    prior = imfilter(prior,kern);
    prior = prior ./ trapz(prior(~isnan(prior)));

    
    % Likelihood
    likelihood = prod(pf_power,3) .* exp(-params.timeBin .* ...
                                         sum(pf_map,3));
    likelihood = likelihood ./ trapz(likelihood(~ ...
                                                isnan(likelihood)));
    
    % Calculate posteriors
    post_u = prior .* likelihood;
    post_u(isnan(post_u)) = 0;
    post_u = imfilter(post_u,gKern);
    recon.update.posterior(:,:,t) = post_u;
    recon.update.loc(t,:) = getMax(post_u);
    
    post_cc = p_map .* CC .* likelihood;
    post_cc(isnan(post_cc)) = 0;
    post_cc = imfilter(post_cc,gKern);
    recon.cc.posterior(:,:,t) = post_cc;
    recon.cc.loc(t,:) = getMax(post_cc);

    
    post_nocc = p_map .* likelihood;
    post_nocc(isnan(post_nocc)) = 0;
    post_nocc = imfilter(post_nocc,gKern);
    recon.nocc.posterior(:,:,t) = post_nocc;
    recon.nocc.loc(t,:) = getMax(post_nocc);
    
        
    
    %% POP VECTOR
    warning off

    for ii = 1:length(bins) - 1
        for jj = 1:length(bins) - 1
            dist_map(ii,jj) = pdist([n(:)'; squeeze(spike_map(ii,jj,:))'],'cosine');
        end
    end
    recon.pv.loc(t,:) = getMax(1 - dist_map);
    recon.pv.dist_map(:,:,t) = 1 - dist_map;
        
    if live_plot
        cla
        % No CC
        subplot(4,4,1);
        imagesc(likelihood);
        ylabel('Bayesian_noCC');
        title('Likelihood');
        subplot(4,4,2)
        imagesc(CC);
        title('CC')
        subplot(4,4,3)
        title('Posterior');
        imagesc(recon.nocc.posterior(:,:,t));
        ax = subplot(4,4,4);
        hold on
        plot(recon.nocc.loc(t,2)/2,...
             recon.nocc.loc(t,1)/2,'rx');
        plot(loc(t,2),loc(t,1),'bo');
        hold off
        xlim([lower upper])
        ylim([lower upper])
        set(ax,'YDir','reverse');
        
        % CC
        subplot(4,4,5);
        imagesc(likelihood);
        ylabel('Bayesian_noCC');
        title('Likelihood');
        subplot(4,4,6)
        imagesc(p_map);
        title('Prior (p_map)')
        subplot(4,4,7)
        title('Posterior');
        imagesc(recon.cc.posterior(:,:,t));
        ax = subplot(4,4,8);
        hold on
        plot(recon.cc.loc(t,2)/2,...
             recon.cc.loc(t,1)/2,'rx');
        plot(loc(t,2),loc(t,1),'bo');
        hold off
        xlim([lower upper])
        ylim([lower upper])
        set(ax,'YDir','reverse');
        
        % Update
        subplot(4,4,9);
        imagesc(likelihood);
        ylabel('Bayesian_noCC');
        title('Likelihood');
        subplot(4,4,10)
        imagesc(prior);
        title('Prior')
        subplot(4,4,11)
        title('Posterior');
        imagesc(recon.cc.posterior(:,:,t));
        ax = subplot(4,4,12);
        hold on
        plot(recon.update.loc(t,2)/2,...
             recon.update.loc(t,1)/2,'rx');
        plot(loc(t,2),loc(t,1),'bo');
        hold off
        xlim([lower upper])
        ylim([lower upper])
        set(ax,'YDir','reverse');
        
        
        % Population vector
        subplot(4,4,13);
        ylabel('Population Vector');
        subplot(4,4,14)
        imagesc(mean(spike_map,3));
        title('Spike Map')
        subplot(4,4,15)
        title('Distance Map');
        imagesc(1-dist_map);
        ax = subplot(4,4,16);
        hold on
        plot(recon.pv.loc(t,2)/2,...
             recon.pv.loc(t,1)/2,'rx');
        plot(loc(t,2),loc(t,1),'bo');
        hold off
        xlim([lower upper])
        ylim([lower upper])
        set(ax,'YDir','reverse');
        keyboard
    end
end

save('reconstruction.mat','loc','recon','bins');



























































%% SUBFUNCS
function [max_idx] = getMax(map)
[~,argmax_idx] = max(map(:));
[max_idx(1),max_idx(2)] = ind2sub(size(map),argmax_idx);


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



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    