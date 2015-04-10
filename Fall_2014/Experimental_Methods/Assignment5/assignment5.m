function assignment5

%% STEP 1 A/B
close all
rng('shuffle', 'twister');

% Parameters:
dotCorr = .5;
dotSD = .5;
nTrials = 1000;

% Criterion values
rightCrit = linspace(-1,1,100);

% Plot ROCs
% First varying dotSD
axis square
fig(1) = figure(1);
subplot(1,2,1)
additive = [0 .125 .25 .375 .5];
colors = [linspace(1,.5,length(additive))' zeros(length(additive),1) zeros(length(additive),1)];
hold on
for j = 1:length(additive)
        for i = 1:length(rightCrit)
            simHR(i) = sum(normrnd(dotCorr, dotSD+additive(j), 1, nTrials) ...
                           > rightCrit(i)) / nTrials;
            simFA(i) = sum(normrnd(-dotCorr, dotSD+additive(j), 1, nTrials) ...
                           > rightCrit(i)) / nTrials;
        end

        % Plot Simulation
        plot(simFA,simHR, '.', 'Color', colors(j,:), 'MarkerSize', 10);
        
        % Compute and plot theoretical values
        [theoryHR theoryFA] = theoryROC(dotCorr,dotSD+additive(j),rightCrit);
        h(j) = plot(theoryFA, theoryHR, 'Color', colors(j,:), 'LineWidth', 2);
end
title('Variable SD (M = .5)');
xlabel('FA');
ylabel('HR');
legend(h, '0.50 (dotSD)', '0.73', '0.75', '0.88', '1.00', 'Location', 'SE');
set(gca,'FontSize',16);
hold off

% Then varying dotCorr
colors = [zeros(length(additive),1) linspace(1,.5,length(additive))' zeros(length(additive),1)];
subplot(1,2,2)
hold on
for j = 1:length(additive)
    simHR = zeros(1,length(rightCrit));
    simFA = zeros(1,length(rightCrit));
    for i = 1:length(rightCrit)
        % Simulated values
        simHR(i) = sum(normrnd(dotCorr+additive(j), dotSD, 1, nTrials) ...
                       > rightCrit(i)) / nTrials;
        simFA(i) = sum(normrnd(-abs(dotCorr+additive(j)), dotSD, 1, nTrials) ...
                       > rightCrit(i)) / nTrials;
        
    end

    % Plot Simulation
    plot(simFA,simHR, '.', 'Color', colors(j,:), 'MarkerSize', 10);
        
    % Compute and plot theoretical values
    clear theoryHR theoryFA
    [theoryHR theoryFA] = theoryROC(dotCorr+additive(j),dotSD,rightCrit);
    h(j) = plot(theoryFA, theoryHR, 'Color', colors(j,:), 'LineWidth', 2);
end
title('Variable Mean (SD = .5)');
xlabel('FA');
ylabel('HR');
legend(h, '0.50 (dotCorr)', '0.73', '0.75', '0.88', '1.00', 'Location', 'SE');
set(gca,'FontSize',16);
hold off




%% STEP 2
% Calculate payoff given diff outcome weights and probability of stimulus going right

% Parameters
Vh =     [0 1];
Vfa =    [-1 -1];
Vm =     [-2 -1];
Vcr =    [0 3];
pRight = [.5 .6];

fig(2) = figure(2);
for k = 1:2
    for i = 1:length(rightCrit)
        payoff(k,i) = calcPayoff(dotCorr, dotSD, pRight(k), Vh(k), Vfa(k), Vm(k), Vcr(k), rightCrit(i));
    end
end
hold on
plot(rightCrit, payoff(1,:),'r','LineWidth',2);
plot(rightCrit, payoff(2,:),'g','LineWidth',2);
title('Estimated Payoffs');
xlabel('Criterion Resp. "Right"');
ylabel('Payoff');
legend(sprintf('pr(%1.1f) Vh(%d) Vfa(%d) Vm(%d) Vcr(%d)',pRight(1),Vh(1),Vfa(1),Vm(1),Vcr(1)), ...
       sprintf('pr(%1.1f) Vh(%d) Vfa(%d) Vm(%d) Vcr(%d)',pRight(2),Vh(2),Vfa(2),Vm(2),Vcr(2)), ...
       'Location', 'SW');
set(gca,'FontSize',16);
hold off




%% STEP 3
% Calculate "theoretical" and data dPrime
newDotCorrs = linspace(0,1,100);
newDotSDs = linspace(0,1,100);

fig(3) = figure(3);
% Varying dotCorr
hold on
for i = 1:length(newDotCorrs)
    [dPrime(i) dPrimeData(i)] = calcdPrime(newDotCorrs(i), dotSD, rightCrit);
end
h3(1) = plot(dPrime, dPrimeData,'ro','MarkerSize',10);

% Varying dotSD
for i = 1:length(newDotSDs(i))
    [dPrime(i) dPrimeData(i)] = calcdPrime(dotCorr, newDotSDs(i), rightCrit);
end
h3(2) = plot(dPrime, dPrimeData,'g.','MarkerSize',10);
title('d'' vs d''Data');
xlabel('d''');
ylabel('d''Data');
legend(h3, 'dotCorrs 0 -> 1', 'dotSD 0 -> 1', 'Location', 'SE');
set(gca,'Fontsize',16);
hold off




%% STEP 4
fig(4) = figure(4);

% Parameters
rightCrit = 0;
pRight = .75;
dotSD = .5;
dotCorr = linspace(0,1,100);
nTrials = 1000;
nRightTrials = nTrials * pRight

hold on
% Simulate Y/N Correct Probability
% For each dot correlation
trial = [];
resp = [];
pCorrYN = [];
for i = 1:length(dotCorr)
    trial(i,:) = shuffle([repmat(dotCorr(i),1,nRightTrials) repmat(-dotCorr(i),1,nTrials-nRightTrials)]);
    % For each trial
    for j = 1:length(trial)
        % Simulate response
        resp(i,j) = normrnd(trial(i,j),dotSD) > rightCrit;
    end
    pCorrYN(i,1) = sum(resp(i,:) == (trial(i,:)>0)) / nTrials;
    
    % Calculate ROC
    [HR FA] = theoryROC(dotCorr(i), dotSD, linspace(-1,1,100));
    AREA(i) = trapz(HR,FA)+1;
end
h4(1) = plot(dotCorr,pCorrYN','.r','MarkerSize',10);

% Simulate 2AFC Correct Probability
respSig = [];
respNoise = [];
pCorrAFC = [];
% For each dot correlation
for i = 1:length(dotCorr)
    % Simulate a response to each dot patch
    respSig(i,:) = normrnd(dotCorr(i),dotSD,1,nTrials);   % "Signal Patch"
    respNoise(i,:) = normrnd(-dotCorr(i),dotSD,1,nTrials);  % "Noise Patch"
    pCorrAFC(i,1) = sum(respSig(i,:) > respNoise(i,:)) / nTrials;
    
    % Calculate d'
    dPrimeAFC(i) = abs((dotCorr(i) - (-dotCorr(i)))/dotSD);
end
h4(2) = plot(dotCorr,pCorrAFC','.g','MarkerSize',10);
h4(3) = plot(dotCorr,AREA,'g','LineWidth',2);
title('p(Correct) in 2AFC vs Y/N Simulations');
xlabel('Dot Correlations');
ylabel('p(Correct)');
legend(h4, 'Y/N', '2AFC','ROC Area', 'Location', 'SE');
set(gca,'FontSize',16);
hold off

% Plot d' vs Percent Correct
fig(5) = figure(5);
plot(dPrimeAFC,pCorrAFC','LineWidth',2);
title('2AFC d'' vs Percent Correct');
xlabel('d''');
ylabel('p(Correct)');
set(gca,'FontSize',16);

% Print out figs to file
for i = 1:length(fig)
    print(fig(i), sprintf('A5_Fig%d.pdf',i), '-dpdf');
end


    

return












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [dPrime dPrimeData] = calcdPrime(dotCorr, dotSD, rightCrit)
% Calculate dPrime for a given dotCorr and dotSD

% First get HR and FR
[HR FA] = theoryROC(dotCorr, dotSD, rightCrit);

% Then calculate d' for data
dPrimeData = mode(abs(norminv(HR) - norminv(FA)));

% Calculate theoretical dPrime?
dPrime = abs((dotCorr - (-dotCorr))/dotSD);

return



function payoff = calcPayoff(dotCorr, dotSD, pRight, Vh, Vfa, Vm, Vcr, rightCrit)
% Calculate payoff given payoff for outcomes, probability of a right-stimulus, dotCorr, dotSD

% First get HR and FR
[HR FA] = theoryROC(dotCorr, dotSD, rightCrit);
%keyboard

% Prior odds of stim:
pS = pRight;
pN = 1 - pRight;

% Probabilities of different outcomes (probability of a certain response given conditional prob)
pHit = HR * pS;
pFA = FA * pN;
pMiss = (1 - HR) * pS;
pCR = (1 - FA) * pN;

% Calculate payoff (sum of weighted outcome probabilities)?
payoff =  (pHit * Vh) + (pFA * Vfa) + (pMiss * Vm) + (pCR * Vcr);

return



function [HR FA] = theoryROC(dotCorr, dotSD, rightCrit)
% Calculate theoretical HR and FA
% (normcdf calculates the position in the probability distribution of the criterion)

HR = 1 - normcdf(rightCrit,repmat(dotCorr,1,length(rightCrit)),repmat(dotSD,1,length(rightCrit)));
FA = 1 - normcdf(rightCrit,repmat(-dotCorr,1,length(rightCrit)),repmat(dotSD,1,length(rightCrit)));

return











