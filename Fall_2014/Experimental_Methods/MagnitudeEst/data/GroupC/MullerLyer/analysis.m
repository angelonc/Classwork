%% MULLER-LYER PLOTTING AND ANALYSIS

% Directories
mainDir = '/Users/chrisan/Documents/MagnitudeEst/data/GroupC/MullerLyer/';
dataDir = dir([mainDir 'sub*']);

n = length(dataDir);

% Build data structure
dat = [];
for i = 1:lengthn
    str = dir([mainDir dataDir(i).name '/*.csv']);
    dat{i} = [repmat(i,24,1) csvread([mainDir dataDir(i).name '/' str.name],1)];
end
dat = vertcat(dat{:});

% Look at in vs out arrows across all conds
for i = 1:n
    inMean(i,:) = mean(dat(dat(:,1) == i & dat(:,5) > 90,6:9));
    inSD(i,:) = std(dat(dat(:,1) == i & dat(:,5) > 90,6:9));
    outMean(i,:) = mean(dat(dat(:,1) == i & dat(:,5) < 90,6:9));
    outSD(i,:) = std(dat(dat(:,1) == i & dat(:,5) < 90,6:9));
end
inAllMean = mean(mean(inMean));
inSEM = mean(mean(inSD))/sqrt(n);
outAllMean = mean(mean(outMean));
outSEM = mean(mean(outSD))/sqrt(n);
barwitherr([outAllMean], [outSEM])






    
