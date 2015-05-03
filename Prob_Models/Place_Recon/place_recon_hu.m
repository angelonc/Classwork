function [loc, re_loc, post, re_loc_ncc, post_nocc, pv_loc, l, u] = place_recon_hu(pv);

%% Bayesian Reconstruction of Spatial Location

load delivery_subject1_session1.mat

%% PREPROCESSING

rescale = 36/1.7;

% Extract relevant info
x = [sessionData.behav.x];
y = [sessionData.behav.y];
sr = 1/sessionData.epochSize;

% Calculate velocity
v = abs([0 diff(sqrt(x.^2 + y.^2))]); % vu/s

% Align spikes
spikes = sessionData.fr' ./ rescale;

% Index still periods
move_idx = v > .001;


%% MAPS

% Make bins
u = ceil(max([max(x(:)) max(y(:))]));
l = floor(min([min(x(:)) min(y(:))]));
nBins = 70;
bins = linspace(l,u,nBins + 1);
[X1,Y1] = meshgrid(bins,bins);
edges{1} = bins;
edges{2} = bins;
params.binSize = mean(diff(bins));

% Make gaussian kernel
sigma = (12 / 2*sqrt(2*log(2)))/params.binSize;
G_kern = fspecial('gaussian',[7*ceil(sigma) 7*ceil(sigma)],sigma);


% Index training period (first 75%)
train_idx = zeros(1,length(x));
train_idx(1:round(length(x)*.75)) = 1;
train_idx = logical(train_idx);

v_train = v(train_idx & move_idx);
x_train = x(train_idx & move_idx);
y_train = y(train_idx & move_idx);
spikes_train = spikes(:,train_idx & move_idx);

% Make sampling and sampling probability map
[N,C] = hist3([x_train' y_train'],'Edges',edges);
s_map = N(1:end-1,1:end-1) / sr; % Convert to seconds
p_map = s_map / sum(s_map(:));
s_map = imfilter(s_map,G_kern,'same');
s_map(s_map == 0) = NaN;

disp('Activity maps...');
f_map = zeros(70,70,32);
for i = 1:size(spikes,1)
    % Make activity map for each cell
    f_map(:,:,i) = genFiringFields(x_train,y_train,spikes_train(i,: ...
                                                      ),bins);
    f_map_s(:,:,i) = imfilter(f_map(:,:,i),G_kern);
    
    % Smooth and divide by sampling map to make place field map
    pf_map(:,:,i) = f_map_s(:,:,i) ./ s_map;
end


%% RECONSTRUCTIONS
disp('Reconstruction...');
% Set up time bins (last 25%)

params.timeBin = 1;
samps_per_bin = round(sr * params.timeBin);
n_bins = floor(length(x(~train_idx))/4);

x_test = x(~train_idx); 
y_test = y(~train_idx);
v_test = v(~train_idx);
spikes_test = spikes(:,~train_idx);

%f_map = zeros(70,70,size(spikes,1));
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
        keyboard

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
    if pv
        for ii = 1:length(bins)-1
            for jj = 1:length(bins) - 1
                dist_map(ii,jj) = pdist([n(:)'; squeeze(f_map_s(ii,jj,:))'],'cosine');
            end
        end
        [~,argmax_idx] = min(dist_map(:));
        [pv_loc(t,1),pv_loc(t,2)] = ind2sub(size(dist_map),argmax_idx);
    end
    

%      ax = subplot(1,4,1);
%      plot(loc(t,2),loc(t,1),'o');
%      xlim([l u])
%      ylim([l u])
%      set(ax,'YDir','reverse');
%      subplot(1,4,2)
%      imagesc(CC)
%      subplot(1,4,3)
%      imagesc(post{t})
%      subplot(1,4,4)
%      imagesc(1 - dist_map)
%      keyboard
    
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
    
            % Make an activity map for time interval
        %if pv
%    f_map(:,:,unit) = genFiringFields(x_test(t_bin(t,:)),y_test(t_bin(t,:)), ...
%                                             spikes_test(unit,t_bin(t,:)),bins);
%       end
        %f_map(:,:,unit) = imfilter(f_map(:,:,unit),G_kern);
end
