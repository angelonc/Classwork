function assignment1

% Experimental Methods Assignment 1

clear all
close all

% Data and conditions
if ~exist('data','var') data = importdata('~/Documents/Classwork/Fall_2014/Experimental_Methods/Assignment1/MullerRawData2014.txt'); end
f1 = logical([1 0 1 1 0 0 1 0]); % line length; 1 = 120, 0 = 75
f2 = [0 1 1 2 0 2 0 1]; % arrows, 0 = in, 1 = out, 2 = no arrows
n = length(data);

% Specify a cutoff to remove repeat conditions
cutoff = 6;

% Find difference between apparent length and actual, if specified:
Diff = 1;
dat = nan([size(data)]);
if Diff
    dat(:,f1) = data(:,f1) - 120;
    dat(:,~f1) = data(:,~f1) - 75;
else
    dat = data;
end


% t-tests w. Bonferroni correction
for i = 1:cutoff
    if Diff
        in = 0;
    elseif f1(i)
        in = 120;
    elseif ~f1(i)
        in = 75;
    end
    [h(i), p(i)] = ttest(dat(:,i),in,'Alpha',.05/cutoff);
end

% Summary statistics
means = mean(dat(:,1:cutoff));
SEM = std(dat(:,1:cutoff))/sqrt(n);

% Format data to plot for bar function
pMeans = zeros(length(unique(f1)),length(unique(f2)));
pSEM = zeros(length(unique(f1)),length(unique(f2)));
for i = 1:size(pMeans,1)
    for j = 1:size(pMeans,2)
        pMeans(i,j) = means(f1(1:6)==(i-1) & f2(1:6)==(j-1));
        pSEM(i,j) = SEM(f1(1:6)==(i-1) & f2(1:6)==(j-1));
        pH(i,j) = h(f1(1:6)==(i-1) & f2(1:6)==(j-1));
    end
end
for i = 1:2
    for j = 1:3
        if pH(i,j) in = '*'; else in = ' '; end
        hString{i,j} = in;
    end
end


% Short line
close all
ax(1) = subplot(1,2,1);
hold on
hb = bar(ax(1), pMeans(1,:),'g');
he = errorbar(ax(1), pMeans(1,:),pSEM(1,:),'.k');
set(gca, 'XTick', 1:cutoff/2, 'XTickLabel', {'>-<', '<->', '-'});
title('75 Line -- Muller-Lyer Illusion');
xlabel('Stimulus Type')
ylabel('Apparent Length - Actual')
if ~Diff
    plot([0 max(hb.XData)+1], [75 75],'--k')
    text(.1, 76, '(Actual Line Length)')
end
% Notate significant bars
text(he.XData-.075,he.YData+pSEM(1,:),hString(1,:), 'FontSize', 16);

% Long line
ax(2) = subplot(1,2,2);
hold on
hb = bar(ax(2), pMeans(2,:),'r');
he = errorbar(ax(2), pMeans(2,:),pSEM(2,:),'.k');
set(gca, 'XTick', 1:cutoff/2, 'XTickLabel', {'>--<', '<-->', '--'});
title('120 Line -- Muller-Lyer Illusion');
xlabel('Stimulus Type')
ylabel('Apparent Length - Actual')
if ~Diff
    plot([0 max(hb.XData)+1], [120 120],'--k')
    text(.1, 122, '(Actual Line Length)')
end
text(he.XData-.075,he.YData+pSEM(2,:),hString(2,:), 'FontSize', 16);
hold off
% Scale both subplots to the larger axis limits
linkaxes([ax(2) ax(1)], 'xy')

saveas(gcf, ['~/Documents/Classwork/Fall_2014/Experimental_Methods/' ...
             'Assignment1/assignment1.pdf'], 'pdf')

return

