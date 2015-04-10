function [threshSub coherences] = analyzeDots(subs)
% Analyzes data from Audio/Visual Dot Motion Experiment
%   s = subject numbers
subs = [6 7];
for s = 1:length(subs)
    subStruct = importdata(strcat('../data2analyze/sub_',num2str(subs(s)),'.dat'));
    %subStruct = importdata(strcat('../data2analyze/sub_',num2str(subs),'.dat'));
    subData=sortrows(subStruct.data,5);
    subData(subData(:,6)==-1,6)=0; % convert -1's in response column to 0's, to make it easier to calculate probabilities later on
    
    % find unique coherence values
    coherences=unique(subData(:,4));
    
    % split data by sound direction
    organizedData(:,:,1)=subData(1:length(subData)/3,:); % downSound in first sheet (first third of data)
    organizedData(:,:,2)=subData(((length(subData)/3)+1):((2/3)*length(subData)),:); % flatSound in second sheet (second third of data)
    organizedData(:,:,3)=subData((((2/3)*length(subData))+1):length(subData),:); % upSound in third sheet (last third of data)
  
    thresholdData = zeros(length(coherences),3); % where I'll put the probabilities for each coherence, for each sound direction
    
    for d = 1:3 % for all sound directions...
        for c = 1:length(coherences) % for all coherences...
            thisCoherenceIndex = organizedData(:,4,d) == coherences(c); % get rows for this coherence value
            thisCoherence_UPprobability = mean(organizedData(thisCoherenceIndex,6,d)); % get probability of "up" response for this coherence value
            thresholdData(c,d)=thisCoherence_UPprobability;
        end
    end
    
    coherences
    thresholdData
    
    
%     % Plot Threshold Curves (dots)
%     Fig = figure('color',[1 1 1]);
%     plot(coherences,thresholdData(:,1),'.-r','MarkerSize',20);
%     hold
%     plot(coherences,thresholdData(:,2),'.-m','MarkerSize',20);
%     plot(coherences,thresholdData(:,3),'.-b','MarkerSize',20);
%     xlabel('Dot Coherence','FontSize',15);
%     ylabel('Probability of "Up" Response','FontSize',15);
%     legend ('down','flat','up','Location','SouthEast')
%     legend('boxoff')
%     saveas(Fig,strcat('../data2analyze/DotMotion_sub',num2str(subs)),'pdf')
    
%     % Plot Threshold Curves (lines)
%     Fig = figure('color',[1 1 1]);
%     plot(coherences,thresholdData(:,1),'-r','MarkerSize',20);
%     hold
%     plot(coherences,thresholdData(:,2),'-m','MarkerSize',20);
%     plot(coherences,thresholdData(:,3),'-b','MarkerSize',20);
%     xlabel('Dot Coherence','FontSize',15);
%     ylabel('Probability of "Up" Response','FontSize',15);
%     legend ('down','flat','up','Location','SouthEast')
%     legend('boxoff')
    
    % Concatenate
    threshSub(:,:,s) = thresholdData;
end

% Plot subjects
h = figure;
hold on
set(gca, 'FontSize', 15)
errorbar(coherences,mean(mean(threshSub,2),3), std(mean(threshSub,2),0,3)/sqrt(size(threshSub,3)), '--', 'Color', [.5 .5 .5]);
errorbar(coherences,mean(threshSub(:,1,:),3),std(threshSub(:,1,:),0,3)/sqrt(size(threshSub,3)),'.-b','MarkerSize',20);
errorbar(coherences,mean(threshSub(:,2,:),3),std(threshSub(:,2,:),0,3)/sqrt(size(threshSub,3)),'.-m','MarkerSize',20);
errorbar(coherences,mean(threshSub(:,3,:),3),std(threshSub(:,3,:),0,3)/sqrt(size(threshSub,3)),'.-r','MarkerSize',20);
title('Impact of Sound Direction on Motion Perception');
xlabel('Dot Coherence');
ylabel('P(Respond Up)');
xlim([-1 1]);
legend ('Average','Decreasing Tone','Flat Tone','Increasing Tone','Location','SouthEast')
legend('boxoff')
hold off


% Plot subjects - no error
h = figure('color',[1 1 1]);
hold on
set(gca, 'FontSize', 15)
plot(coherences,mean(mean(threshSub,2),3), '--', 'Color', [.5 .5 .5]);
plot(coherences,mean(threshSub(:,1,:),3),'.-b','MarkerSize',20);
plot(coherences,mean(threshSub(:,2,:),3),'.-m','MarkerSize',20);
plot(coherences,mean(threshSub(:,3,:),3),'.-r','MarkerSize',20);
title('Impact of Sound Direction on Motion Perception');
xlabel('Dot Coherence');
ylabel('P(Respond Up)');
xlim([-1 1]);
legend ('Average','Decreasing Tone','Flat Tone','Increasing Tone','Location','SouthEast')
legend('boxoff')
hold off



