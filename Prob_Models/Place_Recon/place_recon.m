function place_recon;

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
spike_idx = alignSpikeTimes(ts,unit);

% Exclude still periods
move_idx = abs(v) > params.minV;
x = x(move_idx);
y = y(move_idx);
ts = ts(move_idx);
spikes = spike_idx(:,move_idx);


%% MAPS

% Make bins
bins = 0:params.binSize:35;
edges{1} = bins;
edges{2} = bins;

% Make gaussian kernel
sigma = params.fwhm/params.binSize;
G_kern = fspecial('gaussian',[3*sigma 3*sigma],sigma);

% Make sampling map
[N,C] = hist3([x' y'],'Edges',edges);
s_map = N(1:end-1,1:end-1) / sr; % Convert to seconds
s_map = imfilter(s_map,G_kern);

% Make activity map for each cell
for i = 1:size(spikes,1)
    f_map{i} = genFiringFields(x,y,spikes(i,:),bins);
    

    
    % Smooth and divide by sampling map to make place field map
    p_map{i} = imfilter(f_map{i},G_kern) ./ s_map;
    imagesc(p_map{i});
        keyboard;
end





    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
        