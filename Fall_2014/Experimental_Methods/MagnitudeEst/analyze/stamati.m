%% MULLER-LYER PLOTTING AND ANALYSIS

%% Retrieve Data
% Directories
mainDir = '~/Documents/Classwork/Experimental_Methods/MagnitudeEst/data/GroupC/MullerLyer/';
dataDir = dir([mainDir 'sub*']);

n = length(dataDir);

% Build data structure
dat = [];
for i = 1:n
    str = dir([mainDir dataDir(i).name '/*.csv']);
    dat{i} = [repmat(i,24,1) csvread([mainDir dataDir(i).name '/' str.name],1)];
end
% Convert to one big matrix
dat = vertcat(dat{:});

%% Analyse Data
% Classic effect
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
%barwitherr([outAllMean], [outSEM])



%% Look at outward arrows of same total length
% total lenght of 140
for i = 1:n
    meanAlpha_short(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 & dat(:,2) == 100, 6:9));
    meanTheta_short(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 & dat(:,2) == 100, 6:9));
    SDAlpha_short(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 & dat(:,2) == 100, 6:9));
    SDTheta_short(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 & dat(:,2) == 100, 6:9));
end
meanAllAlpha_short = mean(mean(meanAlpha_short));
SEMAlpha_short = mean(mean(SDAlpha_short))/sqrt(n);
meanAllTheta_short = mean(mean(meanTheta_short));
SEMTheta_short = mean(mean(SDTheta_short))/sqrt(n);
    
% total length of 180
for i = 1:n
    meanAlpha_mid(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 & dat(:,2) == 140, 6:9));
    meanTheta_mid(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140, 6:9));
    SDAlpha_mid(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 & dat(:,2) == 140, 6:9));
    SDTheta_mid(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140, 6:9));
end
meanAllAlpha_mid = mean(mean(meanAlpha_mid));
SEMAlpha_mid = mean(mean(SDAlpha_mid))/sqrt(n);
meanAllTheta_mid = mean(mean(meanTheta_mid));
SEMTheta_mid = mean(mean(SDTheta_mid))/sqrt(n);

% total length of 220
for i = 1:n
    meanAlpha_long(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 & dat(:,2) == 180, 6:9));
    meanTheta_long(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 & dat(:,2) == 180, 6:9));
    SDAlpha_long(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 & dat(:,2) == 180, 6:9));
    SDTheta_long(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 & dat(:,2) == 180, 6:9));
end
meanAllAlpha_long = mean(mean(meanAlpha_long));
SEMAlpha_long = mean(mean(SDAlpha_long))/sqrt(n);
meanAllTheta_long = mean(mean(meanTheta_long));
SEMTheta_long = mean(mean(SDTheta_long))/sqrt(n);

% Plot
figure;
errorbar(1,meanAllAlpha_short,SEMAlpha_short,'o');
axis([0 7 0 200])
hold all;
errorbar(2,meanAllTheta_short,SEMTheta_short,'o');
hold all;
errorbar(3,meanAllAlpha_mid,SEMAlpha_mid,'o');
hold all;
errorbar(4,meanAllTheta_mid,SEMTheta_mid,'o');
hold all;
errorbar(5,meanAllAlpha_long,SEMAlpha_long,'o');
hold all;
errorbar(6,meanAllTheta_long,SEMTheta_long,'o');

%% Look at the inward facing arrows of the same length and see if angle influenced subject perception of line magnitudes

% total lenght of 100
for i = 1:n
    meanAlpha2_short(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 40 & dat(:,2) == 100, 6:9));
    meanAlpha3_short(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 23.094010000000000 & dat(:,2) == 100, 6:9));
    
    meanTheta2_short(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 23.094010000000000 & dat(:,2) == 100, 6:9));
    meanTheta3_short(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 40 & dat(:,2) == 100, 6:9));
    
    SDAlpha2_short(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 40 & dat(:,2) == 100, 6:9));
    SDAlpha3_short(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 23.094010000000000 & dat(:,2) == 100, 6:9));
    
    SDTheta2_short(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 40 & dat(:,2) == 100, 6:9));
    SDTheta3_short(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 23.094010000000000 & dat(:,2) == 100, 6:9));
end
meanAllAlpha2_short = mean(mean(meanAlpha2_short));
meanAllAlpha3_short = mean(mean(meanAlpha3_short));

meanAllTheta2_short = mean(mean(meanTheta2_short));
meanAllTheta3_short = mean(mean(meanTheta3_short));

SEMAlpha2_short = mean(mean(SDAlpha2_short))/sqrt(n);
SEMAlpha3_short = mean(mean(SDAlpha3_short))/sqrt(n);

SEMTheta2_short = mean(mean(SDTheta2_short))/sqrt(n);
SEMTheta3_short = mean(mean(SDTheta3_short))/sqrt(n);
    
% total length of 140
for i = 1:n
    meanAlpha2_mid(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 40 & dat(:,2) == 140, 6:9));
    meanAlpha3_mid(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140, 6:9));
    
    meanTheta2_mid(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140, 6:9));
    meanTheta3_mid(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 40 & dat(:,2) == 140, 6:9));
    
    SDAlpha2_mid(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 40 & dat(:,2) == 140, 6:9));
    SDAlpha3_mid(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140, 6:9));
    
    SDTheta2_mid(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 40 & dat(:,2) == 140, 6:9));
    SDTheta3_mid(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140, 6:9));
end
meanAllAlpha2_mid = mean(mean(meanAlpha2_mid));
meanAllAlpha3_mid = mean(mean(meanAlpha3_mid));

meanAllTheta2_mid = mean(mean(meanTheta2_mid));
meanAllTheta3_mid = mean(mean(meanTheta3_mid));

SEMAlpha2_mid = mean(mean(SDAlpha2_mid))/sqrt(n);
SEMAlpha3_mid = mean(mean(SDAlpha3_mid))/sqrt(n);

SEMTheta2_mid = mean(mean(SDTheta2_mid))/sqrt(n);
SEMTheta3_mid = mean(mean(SDTheta3_mid))/sqrt(n);

% total length of 180
for i = 1:n
    meanAlpha2_long(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 40 & dat(:,2) == 180, 6:9));
    meanAlpha3_long(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 23.094010000000000 & dat(:,2) == 180, 6:9));
    
    meanTheta2_long(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 23.094010000000000 & dat(:,2) == 180, 6:9));
    meanTheta3_long(i,:) = mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 40 & dat(:,2) == 140, 6:9));
    
    SDAlpha2_long(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 40 & dat(:,2) == 180, 6:9));
    SDAlpha3_long(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,4) == 23.094010000000000 & dat(:,2) == 180, 6:9));
    
    SDTheta2_long(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 40 & dat(:,2) == 180, 6:9));
    SDTheta3_long(i,:) = std(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,4) == 23.094010000000000 & dat(:,2) == 180, 6:9));
end
meanAllAlpha2_long = mean(mean(meanAlpha2_long));
meanAllAlpha3_long = mean(mean(meanAlpha3_long));

meanAllTheta2_long = mean(mean(meanTheta2_long));
meanAllTheta3_long = mean(mean(meanTheta3_long));

SEMAlpha2_long = mean(mean(SDAlpha2_long))/sqrt(n);
SEMAlpha3_long = mean(mean(SDAlpha3_long))/sqrt(n);

SEMTheta2_long = mean(mean(SDTheta2_long))/sqrt(n);
SEMTheta3_long = mean(mean(SDTheta3_long))/sqrt(n);

% Plot
figure;
errorbar(1,meanAllAlpha2_long,SEMAlpha2_long,'o');
axis([0 5 0 200])
hold all;
errorbar(2,meanAllTheta2_long,SEMTheta2_long,'o');
hold all;
errorbar(3,meanAllAlpha3_long,SEMAlpha3_long,'o');
hold all;
errorbar(4,meanAllTheta3_long,SEMTheta3_long,'o');


%% Look at total lenghts for all types of arrows

% total length of 140
mean_out_short = mean([meanAllAlpha_short;meanAllTheta_short]);

for i = 1:n
    SD_out_short(i,:) = std(dat(dat(:,1) == i & (dat(:,5) == 120 & dat(:,4) == 40 & dat(:,2) == 100) ...
        | (dat(:,5) == 150 & dat(:,4) == 23.094010000000000 & dat(:,2) == 100), 6:9));
end

SEM_out_short = mean(mean(SD_out_short))/sqrt(n);


mean_in_mid = mean([meanAllAlpha2_mid;meanAllAlpha3_mid; ...
    meanAllTheta2_mid;meanAllTheta3_mid]);

for i = 1:n
    SD_in_mid(i,:) = std(dat(dat(:,1) == i & (dat(:,5) == 30 & dat(:,4) == 40 & dat(:,2) == 140) ...
        | (dat(:,5) == 30 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140) ...
        | (dat(:,5) == 60 & dat(:,4) == 40 & dat(:,2) == 140) ...
        | (dat(:,5) == 60 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140), 6:9));
end

SEM_in_mid = mean(mean(SD_in_mid))/sqrt(n);


% total length of 180

mean_out_mid = mean([meanAllAlpha_mid;meanAllTheta_mid]);

for i = 1:n
    SD_out_mid(i,:) = std(dat(dat(:,1) == i & (dat(:,5) == 120 & dat(:,4) == 40 & dat(:,2) == 140) ...
        | (dat(:,5) == 150 & dat(:,4) == 23.094010000000000 & dat(:,2) == 140), 6:9));
end

SEM_out_mid = mean(mean(SD_out_mid))/sqrt(n);


mean_in_long = mean([meanAllAlpha2_long;meanAllAlpha3_long; ...
    meanAllTheta2_long;meanAllTheta3_long]);

for i = 1:n
    SD_in_long(i,:) = std(dat(dat(:,1) == i & (dat(:,5) == 30 & dat(:,4) == 40 & dat(:,2) == 180) ...
        | (dat(:,5) == 30 & dat(:,4) == 23.094010000000000 & dat(:,2) == 180) ...
        | (dat(:,5) == 60 & dat(:,4) == 40 & dat(:,2) == 180) ...
        | (dat(:,5) == 60 & dat(:,4) == 23.094010000000000 & dat(:,2) == 180), 6:9));
end

SEM_in_long = mean(mean(SD_in_long))/sqrt(n);

% Plot
figure;
errorbar(1,mean_out_short,SEM_out_short,'o');
axis([0 5 0 200])
hold all;
errorbar(2,mean_in_mid,SEM_in_mid,'o');
hold all;
errorbar(3,mean_out_mid,SEM_out_mid,'o');
hold all;
errorbar(4,mean_in_long,SEM_in_long,'o');