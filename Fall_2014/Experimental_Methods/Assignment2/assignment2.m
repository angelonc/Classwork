function assignment2

%% Assignment 2
% Task: determine the precision to which the slope of a line can be determined using bootstrapping
clear all; close all;
rng('shuffle');

%% Parameters
% True data params
x = [10 20 30 40 50];
trueSlope = 3;
y = trueSlope*x;

% "Observed" data params
noiseSD = [10 20];
nReps = [10 25 50];

% Simulation params
nSimulations = 1000;
nBootstraps = 100;

params = v2struct(x, y, trueSlope, noiseSD, nReps, nSimulations, nBootstraps);

%% For each variation in data parameters
loopCnt = 0;
f = figure;
set(f, 'Position', [0 0 1200 800])
for m = 1:length(noiseSD)
    for l = 1:length(nReps)
        bootSlope = [];
        bootCI = [];
        obsSlope = [];
        theoryCI = [];
        for k = 1:nSimulations
            % Create some data
            cnt = 0;
            for j = 1:length(x)
                for i = 1:nReps(l)
                    cnt = cnt + 1;
                    xObs(cnt) = x(j);
                    yObs(cnt) = y(j) + normrnd(0, noiseSD(m));
                end
            end
            % Line fit of original data
            [B, BINT] = regress(yObs', xObs', .05);
            obsSlope(k) = B(1);
            theoryCI(k,:) = BINT(1,:);
            
            %% Bootstrapped line fits:
            for j = 1:nBootstraps
                % Draw a new sample with replacement
                ind = unidrnd(cnt,1,cnt);
                % Fitting a line, constrained through origin
                bootSlope(k,j) = regress(yObs(ind)',xObs(ind)');
            end
            % Bootstrapped CI
            bootCI(k,:) = prctile(bootSlope(k,:),[2.5 97.5]);

            % Check fit w. plot
            if 1==2
                hold on
                scatter(xObs,yObs)
                yHat = obsSlope(k).*xObs;
                plot([0 xObs],[0 yat])
                xlim([0 60])
                hold off
            end
        end
        %% Plotting
        loopCnt = loopCnt + 1;
        subplot(length(noiseSD), length(nReps), loopCnt)
        hold on
        
        % Plot Data
        [n,bins] = hist(mean(bootSlope,2));
        bar(bins,n,1)
        bootCI_Mean = mean(bootCI);
        theoryCI_Mean = mean(theoryCI);
        for p = 1:2
            h(1) = plot(repmat(bootCI_Mean(p),1,2), [0 max(ylim)], 'r','LineWidth',3);
            h(2) = plot(repmat(theoryCI_Mean(p),1,2), [0 max(ylim)], 'g','LineWidth',3);
        end
        h(3) = plot(repmat(trueSlope,1,2), [0 max(ylim)], ':c', 'LineWidth', ...
                    3);
        
        % Formatting
        %set(gca,'FontSize',14)
        title(sprintf('CI Estimation: Noise SD %02d, N = %02d\n(%d Simulations)', ...
                      noiseSD(m), nReps(l), nSimulations))
        xlabel(sprintf('Sample-Mean Slope Estimates\n(%d Bootstraps)', ...
                       nBootstraps));
        ylabel('Frequency');
        %xlim([floor(min([bins bootCI_Mean theoryCI_Mean])*10)/10 ...
%      ceil(max([bins bootCI_Mean theoryCI_Mean])*10)/10]);
        xlim([2.7 3.3])
        hold off
    end
end

hold on
% Add single legend
hl = legend(h,'Bootstrapped CI', 'Theoretical CI', 'True Mean');
set(hl, 'Position', [.01 .01 .12 .06], 'Units', 'normalized');
hold off

print -painters -dpdf -r300 test.pdf


return


