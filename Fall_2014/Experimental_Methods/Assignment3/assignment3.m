function assignment3

%% FREQUENCY OF SEEING (Pt1)
% Predict a frequency of seeing curve based on a simple model
% Assumptions of the model:
%  - yes = absorptions>crit no = absorptions<crit
%  - Num photons is Poisson distributed
clear all; close all

%% Model setup:
% Data
stimIntensity = [0:3:21];

% Parameters
detectThresh = 9;
nTrials = 10;


%% Simulated FoS Curve
% For each stimulus intensity:
for i = 1:length(stimIntensity)
    % For each trial:
    for j = 1:nTrials
        % Does rand Poiss. sample at intensity level exceed detection threshold?
        simResponse(i,j) = poissrnd(stimIntensity(i)) > ...
            detectThresh;
    end
end


%% Theoretical FoS Curve
% Do a fine-grain simulation of the curve:
stimLevel = [0:.01:21];
for i = 1:length(stimLevel)
    pSeeing(i) = 1 - poisscdf(detectThresh, stimLevel(i));
end


%% Find 75% Threshold
[minVal, ind] = min(abs(pSeeing-.75));


%% Plot
hold on
set(gca,'FontSize',16)
h(1) = plot(stimLevel,pSeeing,'b','LineWidth',3);
h(2) = plot(stimIntensity,mean(simResponse,2),'.r','MarkerSize',20);
h(3) = plot(repmat(stimLevel(ind),1,2), [0 max(ylim)],'--k','LineWidth',2);
title(sprintf('Simulated FoS Curve\n(%02d Trials, Crit = %d)',nTrials, ...
              detectThresh));
xlabel('Stimulus Intensity (# Photons Absorbed)');
ylabel('Predicted p(Seen)');
legend(h,'Theoretical','Simulated','75% p(Seen)','Location', ...
       'East');
hold off
print -painters -dpdf -r300 assignment3.pdf

return



        
        

