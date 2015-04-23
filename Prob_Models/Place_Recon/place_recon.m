function [loc, re_loc, post, pv_loc] = place_recon(pv);

%% Bayesian Reconstruction of Spatial Location

load Neutral_s1_d1.mat

%% PREPROCESSING

% Extract relevant info
x = pos.rawX;
y = pos.rawY;
ts = pos.ts;

% Calculate sample rate of position data
sr = 1000/mean(diff(ts));

% Calculate velocity
v = [0 diff(sqrt(x.^2 + y.^2))*1000]; % cm/s

% Align spikes
spikes = alignSpikeTimes(ts,unit);

% Index still periods
move_idx = abs(v) > params.minV;


%% MAPS

% Make bins
bins = 0:params.binSize:35;
[X1,Y1] = meshgrid(bins,bins);
edges{1} = bins;
edges{2} = bins;

% Make gaussian kernel
sigma = (params.fwhm / 2*sqrt(2*log(2)))/params.binSize;
G_kern = fspecial('gaussian',[7*ceil(sigma) 7*ceil(sigma)],sigma);

% Make circular mask
cx = 35.5; cy = 35.5;
r = 35.5;
ix = 71;
iy = 71;
[mx,my]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
c_mask=((mx.^2+my.^2)<=r^2);
mask = +c_mask(1:end-1,1:end-1);
mask(mask == 0) = nan;

% Index training period (first 75%)
train_idx = zeros(1,length(ts));
train_idx(1:round(length(ts)*.75)) = 1;
train_idx = logical(train_idx);

v_train = v(train_idx & move_idx);
x_train = x(train_idx & move_idx);
y_train = y(train_idx & move_idx);
spikes_train = spikes(:,train_idx & move_idx);

% Make sampling and sampling probability map
[N,C] = hist3([x_train' y_train'],'Edges',edges);
s_map = N(1:end-1,1:end-1) / sr; % Convert to seconds
p_map = s_map / sum(s_map(:)) .* mask;
s_map = imfilter(s_map,G_kern,'same');
s_map = s_map .* mask;

disp('Activity maps...');
for i = 1:size(spikes,1)
    % Make activity map for each cell
    f_map{i} = genFiringFields(x_train,y_train,spikes_train(i,:),bins);
    
    % Smooth and divide by sampling map to make place field map
    pf_map(:,:,i) = imfilter(f_map{i},G_kern) ./ s_map;
    
    % Find place field peak ("center")
    tmp = max(max(pf_map(:,:,i)));
    [pf_max(i,1),pf_max(i,2)] = find(pf_map(:,:,i) == tmp);
end


%% RECONSTRUCTIONS
disp('Reconstruction...');
% Set up time bins (last 25%)

samps_per_bin = round(sr * params.timeBin);
n_bins = floor(length(ts(~train_idx))/30);

x_test = x(~train_idx); 
y_test = y(~train_idx);
v_test = v(~train_idx);
spikes_test = spikes(:,~train_idx);

%f_map = zeros(70,70,size(spikes,1));
n = zeros(1,1,size(spikes,1));
t_bin = zeros(n_bins,samps_per_bin);
f_map = [];
pv_loc = [];
% For each time bin
for t = 1:n_bins
    disp(['Time ' num2str(t)]);
    t_bin(t,:) = [((t-1) * samps_per_bin) + 1 :...
            ((t-1) * samps_per_bin) + samps_per_bin];
    
    % Current location
    loc(t,:) = [mean(x_test(t_bin(t,:)))*2 mean(y_test(t_bin(t,:)))*2];
    
    
    %% BAYESIAN
    % Continuity constraint around previous position
    if t == 1
        CC = ones(70,70);
    else
        % Centered on last position
        %Mu = fliplr([mean(x_test(t_bin(t-1,:))) mean(y_test(t_bin(t-1,:)))]);
        Mu = fliplr(re_loc(t-1,:)/2);
        % Width is present velocity
        Sigma = abs(mean(v_test(t_bin(t,:)))/params.binSize);
        %CC = fspecial('gauss',[ceil(sigma).*7],sigma);
        if Sigma > 60
            Sigma = 60;
        elseif Sigma < 20
            Sigma = 20;
        end
        Sigma = repmat(Sigma,1,2);
        CC = reshape(mvnpdf([X1(:) Y1(:)],Mu,Sigma),71,71);
        CC = CC(1:end-1,1:end-1) ./ sum(CC(:));
    end
   
    % Make an activity map
    for unit = 1:size(spikes,1)
        % Make an activity map for time interval
        if pv
            f_map(:,:,unit) = genFiringFields(x_test(t_bin(t,:)),y_test(t_bin(t,:)), ...
                                              spikes_test(unit,t_bin(t,: ...
                                                              )), ...
                                              bins);
        end
        %f_map(:,:,unit) = imfilter(f_map(:,:,unit),G_kern);
        
        % Get spike count per unit
        n(unit) = sum(spikes_test(unit,t_bin(t,:)));
        pf_power(:,:,unit) = pf_map(:,:,unit) .^ n(unit);
    end
    
    % Calculate posterior
    post{t} = CC .* p_map .* prod(pf_power,3) .* exp(-params.timeBin .* sum(pf_map,3));
    [~,argmax_idx] = max(post{t}(:));
    [re_loc(t,1),re_loc(t,2)] = ind2sub(size(post{t}),argmax_idx);
    
    %% POP VECTOR
    if pv
        for ii = 1:length(bins)-1
            for jj = 1:length(bins) - 1
                clear tmp;
                tmp = corrcoef(squeeze(pf_map(ii,jj,:)),squeeze(f_map(ii,jj,:)));
                corr_map(ii,jj) = tmp(2,1);
            end
        end
        [~,argmax_idx] = max(corr_map(:));
        [pv_loc(t,1),pv_loc(t,2)] = ind2sub(size(corr_map),argmax_idx);
    end
    

% $$$     ax = subplot(1,4,1)
% $$$     p1 = plot(loc(t,2),loc(t,1),'o')
% $$$     xlim([0 70])
% $$$     ylim([0 70])
% $$$     set(ax,'YDir','reverse');
% $$$     subplot(1,4,2)
% $$$     imagesc(CC)
% $$$     subplot(1,4,3)
% $$$     imagesc(post{t})
% $$$     subplot(1,4,4)
% $$$     imagesc(corr_map)
% $$$     keyboard
    
end


    
    

if 0
    % Centered on last position
    Mu = [mean(x_test(t_bin(t-1,:))) mean(y_test(t_bin(t-1,:)))];
    % Width is present velocity
    Sigma = abs(mean(v_test(t_bin(t,:)))/params.binSize);
    %CC = fspecial('gauss',[ceil(sigma).*7],sigma);
    if Sigma > 60
        Sigma = 20/2;
    elseif Sigma < 20
        Sigma = 10/2;
    end
    Sigma = repmat(Sigma,1,2);
    CC = reshape(mvnpdf([X1(:) Y1(:)],Mu,Sigma),71,71);
    CC = CC(1:end-1,1:end-1);
end

    
    
    
    
    
    
    
    














%% SUBFUNCTIONS
    
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

function [pf_idx] = findPlaceFields(p_map,params)
    
% Threshold the place map
fr_thresh = max(p_map(:)) * params.minPFAct;
area_thresh = params.minPFSize / params.binSize;
p_map_thresh = p_map > fr_thresh;

% Find connected regions
CC = bwconncomp(p_map_thresh);

% Extract fields greater than 80 cm-sq
pf_cnt = 0;
for j = 1:length(CC.PixelIdxList)
    
    if length(CC.PixelIdxList{j}) > area_thresh
        pf_cnt = pf_cnt + 1;
        pf_idx{pf_cnt} = CC.PixelIdxList{j};
    end
end
        