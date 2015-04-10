load('~/Documents/Classwork/Statistics_711/ToAskTim.mat')

% Run regression:
[B,~,r] = regress(beh(:,1), [ones(length(lfp(:,1)),1) lfp(:,1)]);

% Check for beh times:
predBeh = [round([ones(length(beh(:,1)),1), beh(:,1)] * B1) beh(:, ...
                                                  1)];


predLFP = (beh(:,a) - B(1) - r) / B(2);
[lfp(:,a) round(predLFP)];
        
alignedLFP = (behTimes(:,a) - B(1)) / B(2);

% Why doesnt the reverse work?
[B1,~,r1] = regress(lfp(:,1), [ones(length(beh(:,1)),1) beh(:,1)]);
[B1,~,r1] = regress(lfp(:,1), beh(:,1));

% Gives error