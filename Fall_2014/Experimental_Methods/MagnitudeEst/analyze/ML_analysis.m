function stats = ML_analysis

%% MULLER-LYER PLOTTING AND ANALYSIS
clear all

% Directories
mainDir = '~/Documents/Classwork/Experimental_Methods/MagnitudeEst/data/GroupC/MullerLyer/';
dataDir = dir([mainDir 'sub*']);
addpath(genpath('~/Documents/MATLAB/anova_rm'));


%% Build data structure
n = length(dataDir);
dat = [];
for i = 1:n
    str = dir([mainDir dataDir(i).name '/*.csv']);
    dat{i} = [repmat(i,24,1) csvread([mainDir dataDir(i).name '/' str.name],1)];
end
dat = vertcat(dat{:});



%% Check basic illusion
% Stats across blocks
for i = 1:n
    inMean(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) > 90,6:9), ...
                          2)',2);
    inSD(i) = std(mean(dat(dat(:,1) == i & dat(:,5) > 90,6:9), ...
                          2));
    outMean(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) < 90,6:9), ...
                           2)',2);
    outSD(i) = std(mean(dat(dat(:,1) == i & dat(:,5) < 90,6:9), ...
                           2));
end

% Summary stats
stats.classic.inM = mean(mean(inMean));
stats.classic.inSEM = std(inMean)/sqrt(n);
stats.classic.outM = mean(mean(outMean));
stats.classic.outSEM = std(outMean)/sqrt(n);
[h,stats.classic.p,ci,stats.classic.stats] = ttest(inMean, outMean);


close all
h(1) = figure(1);
hold on
plotM = [stats.classic.inM stats.classic.outM];
plotSEM = [stats.classic.inSEM stats.classic.outSEM];
plotCol = [1 0 0; 0 0 1];
he = errorbar(plotM,plotSEM,'k.');
for i = 1:length(plotM)
    p(i) = plot(i, plotM(i), '.', 'color', plotCol(i,:), 'markersize', 30);
end
set(gca, 'XTick', 1:2, 'XTickLabel', {'>-<', '<->'});
title(['Classic Muller-Lyer Illusion - Average Estimates per Arrow ' ...
       'Direction']);
xlabel('Arrow Direction')
ylabel('Magnitude Estimates')
legend(p,'>-<','<->', 'location','ne');
set(gca,'FontSize',16)
hold off



%% Look at outward arrows of same total length
% total lenght of 140
for i = 1:n
    meanAlpha_short(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 &...
                                         dat(:,2) == 100, 6:9),2)',2);
    meanTheta_short(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 &...
                                         dat(:,2) == 100, 6:9),2)',2);
    SDAlpha_short(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 &...
                                 dat(:,2) == 100, 6:9),2));
    SDTheta_short(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 &...
                                 dat(:,2) == 100, 6:9),2));
end
meanAllAlpha_short = mean(mean(meanAlpha_short));
SEMAlpha_short = std(meanAlpha_short)/sqrt(n);
meanAllTheta_short = mean(mean(meanTheta_short));
SEMTheta_short = std(meanTheta_short)/sqrt(n);
    
% total length of 180
for i = 1:n
    meanAlpha_mid(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 &...
                                         dat(:,2) == 140, 6:9),2)',2);
    meanTheta_mid(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 &...
                                         dat(:,2) == 140, 6:9),2)',2);
    SDAlpha_mid(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 &...
                                 dat(:,2) == 140, 6:9),2));
    SDTheta_mid(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 &...
                                 dat(:,2) == 140, 6:9),2));
end
meanAllAlpha_mid = mean(mean(meanAlpha_mid));
SEMAlpha_mid =std(meanAlpha_mid)/sqrt(n);
meanAllTheta_mid = mean(mean(meanTheta_mid));
SEMTheta_mid = std(meanTheta_mid)/sqrt(n);

% total length of 220
for i = 1:n
    meanAlpha_long(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 &...
                                         dat(:,2) == 180, 6:9),2)',2);
    meanTheta_long(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 &...
                                         dat(:,2) == 180, 6:9),2)',2);
    SDAlpha_long(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 120 & dat(:,4) == 40 &...
                                 dat(:,2) == 180, 6:9),2));
    SDTheta_long(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 150 & dat(:,4) == 23.094010000000000 &...
                                 dat(:,2) == 180, 6:9),2));
end
meanAllAlpha_long = mean(mean(meanAlpha_long));
SEMAlpha_long = std(meanAlpha_long)/sqrt(n);
meanAllTheta_long = mean(mean(meanTheta_long));
SEMTheta_long = std(meanTheta_long)/sqrt(n);

stats.out.alphaMeans = [meanAllAlpha_short meanAllAlpha_mid ...
                    meanAllAlpha_long];
stats.out.alphaSEM = [SEMAlpha_short SEMAlpha_mid SEMAlpha_long];
stats.out.thetaMeans = [meanAllTheta_short meanAllTheta_mid ...
                    meanAllTheta_long];
stats.out.thetaSEM = [SEMTheta_short SEMTheta_mid SEMTheta_long];
stats.out.ANOVA = anova_rm({[meanAlpha_short' meanAlpha_mid' ...
                    meanAlpha_long'] [meanTheta_short' meanTheta_mid' ...
                    meanTheta_long']},'off');


% Plot
h(2) = figure(2);
hold on
errorbar(stats.out.alphaMeans, stats.out.alphaSEM,'.k');
errorbar(stats.out.thetaMeans, stats.out.thetaSEM,'.k');
p(1) = plot(stats.out.alphaMeans, '.-r', 'markersize', 30);
p(2) = plot(stats.out.thetaMeans, '.-b', 'markersize', 30);
legend(p,'120 Degrees','150 Degrees','location','nw');
title(['Effect of Arrow Angle on >-< Lines Matched for Overall ' ...
       'Length']);
xlabel('Overall Figure Length (pixels)');
ylabel('Estimated Line Length');
text(0,-10,'(error bars == SEM)')
set(gca, 'XTick', 1:4, 'XTickLabel', {'140', '180', '220'});
set(gca,'FontSize',16)
hold off

%% Look at the inward facing arrows of the same length and see if angle influenced subject perception of line magnitudes

% total lenght of 100
for i = 1:n
    meanAlpha2_short(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,2) == 100, 6:9),2)',2);
    meanTheta2_short(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,2) == 100, 6:9),2)',2);
    SDAlpha2_short(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,2) == 100, 6:9),2));
    SDTheta2_short(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,2) == 100, 6:9),2));
end
meanAllAlpha2_short = mean(mean(meanAlpha2_short));
meanAllTheta2_short = mean(mean(meanTheta2_short));

SEMAlpha2_short = std(meanAlpha2_short)/sqrt(n);
SEMTheta2_short = std(meanTheta2_short)/sqrt(n);

    
% total length of 140
for i = 1:n
    meanAlpha2_mid(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,2) == 140, 6:9),2)',2);
    meanTheta2_mid(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,2) == 140, 6:9),2)',2);
    SDAlpha2_mid(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,2) == 140, 6:9),2));
    SDTheta2_mid(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,2) == 140, 6:9),2));
end
meanAllAlpha2_mid = mean(mean(meanAlpha2_mid));
meanAllTheta2_mid = mean(mean(meanTheta2_mid));

SEMAlpha2_mid = std(meanAlpha2_mid)/sqrt(n);
SEMTheta2_mid = std(meanTheta2_mid)/sqrt(n);

% total length of 180
for i = 1:n
    meanAlpha2_long(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,2) == 180, 6:9),2)',2);
    meanTheta2_long(i) = mean(mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,2) == 180, 6:9),2)',2);
    SDAlpha2_long(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 30 & dat(:,2) == 180, 6:9),2));
    SDTheta2_long(i) = std(mean(dat(dat(:,1) == i & dat(:,5) == 60 & dat(:,2) == 180, 6:9),2));
end
meanAllAlpha2_long = mean(mean(meanAlpha2_long));
meanAllTheta2_long = mean(mean(meanTheta2_long));

SEMAlpha2_long = std(meanAlpha2_long)/sqrt(n);
SEMTheta2_long = std(meanTheta2_long)/sqrt(n);

stats.in.alphaMean = [meanAllAlpha2_short meanAllAlpha2_mid ...
                    meanAllAlpha2_long];
stats.in.alphaSEM = [SEMAlpha2_short SEMAlpha2_mid SEMAlpha2_long];
stats.in.thetaMean = [meanAllTheta2_short meanAllTheta2_mid ...
                    meanAllTheta2_long];
stats.in.thetaSEM = [SEMTheta2_short SEMTheta2_mid SEMTheta2_long];
stats.in.ANOVA = anova_rm({[meanAlpha2_short' meanAlpha2_mid' ...
                    meanAlpha2_long'] [meanTheta2_short' meanTheta2_mid' ...
                    meanTheta2_long']},'off');

% Plot
h(3) = figure(3);
hold on
errorbar(stats.in.alphaMean, stats.in.alphaSEM,'.k');
errorbar(stats.in.thetaMean, stats.in.thetaSEM,'.k');
p(1) = plot(stats.in.alphaMean, '.-r', 'markersize', 30);
p(2) = plot(stats.in.thetaMean, '.-b', 'markersize', 30);
legend(p,'30 Degrees','60 Degrees','location','nw');
title(['Effect of Arrow Angle on <-> Lines Matched for Overall ' ...
       'Length']);
xlabel('Overall Figure Length (pixels)');
ylabel('Estimated Line Length');
set(gca, 'XTick', 1:3, 'XTickLabel', {'100', '140', '180'});
set(gca,'FontSize',16)
hold off

%% Look at relationships between matched pairs of in vs outward arrows
for i = 1:n
    % 100 in, 140 out - 140 total
    meanIn100(i) = mean(mean(dat(dat(:,1) == i & dat(:,2) == 100 &...
                             ((dat(:,5) == 120 & dat(:,4) == 40) |...
                              dat(:,5) == 150 & dat(:,4) == ...
                              40/sqrt(3)),6:9),2)',2);
    meanOut140(i) = mean(mean(dat(dat(:,1) == i & dat(:,2) == 140 &...
                              (dat(:,5) == 30 | dat(:,5) == 60),6:9), ...
                              2)',2);
    % 140 in, 180 out - 180 total
    meanIn140(i) = mean(mean(dat(dat(:,1) == i & dat(:,2) == 140 &...
                             ((dat(:,5) == 120 & dat(:,4) == 40) |...
                              dat(:,5) == 150 & dat(:,4) == ...
                              40/sqrt(3)),6:9),2)',2);
    meanOut180(i) = mean(mean(dat(dat(:,1) == i & dat(:,2) == 180 &...
                              (dat(:,5) == 30 | dat(:,5) == 60),6:9),2)',2);
end
stats.pair.inMeans = [mean(mean(meanIn100)) ...
                    mean(mean(meanIn140))];
stats.pair.inSEM = [std(meanIn100)/sqrt(n) std(meanIn140)/sqrt(n)];
stats.pair.outMeans = [mean(mean(meanOut140)) ...
                    mean(mean(meanOut180))];
stats.pair.outSEM = [std(meanOut140)/sqrt(n) std(meanOut180)/ ...
                    sqrt(n)];
stats.pair.ANOVA = anova_rm({[meanIn100' meanOut140'] [meanIn140' ...
                    meanOut180']},'off');
[h, stats.pair.postHoc.p, ci, stats.pair.postHoc.S] = ...
    ttest(mean([meanIn100' meanOut140'],2), mean([meanIn140' meanOut180'],2));


h(4) = figure(4);
hold on
e(1) = errorbar(stats.pair.inMeans, stats.pair.inSEM, '.k');
e(2) = errorbar(stats.pair.outMeans, stats.pair.outSEM, '.k');
p(1) = plot(stats.pair.inMeans, '.-r', 'markersize', 30);
p(2) = plot(stats.pair.outMeans, '.-b', 'markersize', 30);
legend(p,'>-<','<->','location','nw');
title(['Effect of Arrow Direction on Lines Matched for Overall ' ...
       'Length']);
xlabel('Overall Figure Length (pixels)');
ylabel('Estimated Line Length');
text(1,1,'(error bars == SEM)')
set(gca, 'XTick', 1:2, 'XTickLabel', {'140', '180'});
set(gca,'FontSize',16)
hold off
