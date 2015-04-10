function Prob_Modelling_HW2
%% PROBABILISTIC MODELING HW2
close all;

%%%%%%%%%%%%%%%%%%%%%%%%
%% A: Optimal Estimation
%%%%%%%%%%%%%%%%%%%%%%%%

MUx = 10;
SIGx = sqrt(8);
n = 200;
Rx = [0 30];
Ry = [0 20];

%% A.1: Simulate data
% Make sure x is in range
x_sim = -1 * ones(1,n);
while any(x_sim <= Rx(1)) | any(x_sim >= Rx(2))
    x_sim = normrnd(MUx,SIGx,1,n);
end
% Make sure y is in range
y_sim = -1 * ones(1,n);
while any(y_sim <= Ry(1)) | any(y_sim >= Ry(2))
    y_sim = normrnd(x_sim, .5 + .2 * x_sim);
end

fprintf('A.1: Simulation with mean...\n\n');
figure
set(gca,'FontSize',20);
xlim(Rx);
ylim(Ry);
hold on
scatter(x_sim,y_sim,'ob','MarkerFaceColor','b');
h = lsline;
set(h,'LineWidth',2,'Color',[.1 .1 .1],'LineStyle',':');
hold on
title('Simulated Part-Length and Measurement');
xlabel('Part Length x (mm)');
ylabel('Measurement y (mm)');
pause;


%% A.2/3: Derive MAP and LS as a function of y
x = linspace(Rx(1),Rx(2),200);
y = linspace(Ry(1),Ry(2),200);
[iMAP, LS] = calcMAP_LS(x,y,MUx,SIGx);
[iMAP_sim, LS_sim] = calcMAP_LS(x,y_sim,MUx,SIGx);


fprintf('A.2 and A.3: MAP and BLS estimation...\n\n');
figure
set(gca,'FontSize',20);
xlim(Ry);
ylim(Rx);
hold on
h(1) = scatter(y_sim,x_sim,'ob','MarkerFaceColor','b');
h(2) = lsline;
scatter(y_sim,x(iMAP_sim),'g');
scatter(y_sim,LS_sim,'r');
set(h(2),'LineWidth',2,'Color',[.1 .1 .1],'LineStyle',':');
h(3) = plot(y,x(iMAP),'g','LineWidth',1);
h(4) = plot(y,LS,'r','LineWidth',1);
hold off
title('MAP/LS Estimates');
xlabel('Measurement y (mm)');
ylabel('Predicted Part Length x (mm)');
legend(h,'Simulation','Mean','MAP','LS');
pause;
clear all

%% A.4
fprintf(['A.4\nThe estimators do not follow the mean because they ' ...
         'are biased towards the posterior distribution of the data, ' ...
         'which is biased by the prior. The mean, on the other hand, ' ...
         'is predictive only based on the observations over which ' ...
         'it is averaged, while MAP and BLS estimates have some ' ...
         '"knowledge" of the underlying distributions.\n\n']);
pause;

fprintf(['A.5\nThe two estimators differ slightly because they are ' ...
         'estimating based on different facets of the posterior ' ...
         'distribution. MAP estimates rely on local maxima, which ' ...
         'will tend to rely entirely on the skewness of the posterior. ' ...
         'BLS on the other hand, uses the mean of the posterior, ' ...
         'which will give a better estimate of the central tendencies ' ...
         'of the data, causing it to gravitate towards the tail ' ...
         'end of the distribution (as is seen in the graph, the BLS ' ...
         'estimate is slightly shifted towards the "upper" portion ' ...
         'of the data set that exhibits more variance/spread.\n\n']);
pause;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% B - Theory of Signal Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% B.1/2 Compute ROC curve of source in natural environment
lambdaN = 5;
lambdaS = [5 10 20 100];
crit = linspace(0,40,30);
colIdx = [linspace(.25,1,length(lambdaS))' ...
          zeros(1,length(lambdaS))' zeros(1,length(lambdaS))'];

if 0
    figure
    hold on
    x = 0:120;
    plot(poisspdf(x,lambdaN),'r');
    plot(poisspdf(x,lambdaS(3)),'b');
    hold off
    pause
    close
end

fprintf('B.1 and B.2: ROC curves for various Poisson distributed sources...\n\n');
figure
set(gca,'FontSize',20)
hold on
for i = 1:length(lambdaS)
    HR = 1 - poisscdf(crit,repmat(lambdaS(i),1,length(crit)));
    FA = 1 - poisscdf(crit,repmat(lambdaN,1,length(crit)));
    
    h(i) = plot(FA,HR,'Color',colIdx(i,:),'LineWidth',2);
end
hold off
legend(h, '5 rads/s', '10 rads/s', '29 rads/s', '100 rads/s', ...
       'Location','SE');
title('ROC Curves for Different Signals')
xlabel('FA')
ylabel('HR')
ylim([0 1.01])
xlim([-.01 1])
pause

% Fun video!
vid = input('Press 1 to watch video, 0 to skip: ');
fprintf('\n');
if vid
    figure
    clear HR FA
    for j = 1:length(lambdaS)
        crit = linspace(0,2.5*lambdaS(j),50);
        for i = 1:length(crit)
            clear HR FA
            
            subplot(1,2,1)
            x = 1:max(crit);
            xlim([0 max(crit)]);
            hold on
            cla
            plot(poisspdf(x,lambdaN),'r','LineWidth',3);
            plot(poisspdf(x,lambdaS(j)),'b--','LineWidth',3);
            h1 = plot([crit(i) crit(i)], ylim,'g','LineWidth',2);
            title('Distributions and Criterion');
            xlabel(sprintf('Decay Rate (LambdaSignal = %g)',lambdaS(j)));
            ylabel('Probability Density');
            hold off
            
            subplot(1,2,2)
            hold on
            HR(i) = 1 - poisscdf(crit(i),lambdaS(j));
            FA(i) = 1 - poisscdf(crit(i),lambdaN);
            HR_all = 1 - poisscdf(crit,repmat(lambdaS(j),1,length(crit)));
            FA_all = 1 - poisscdf(crit,repmat(lambdaN,1,length(crit)));
            plot(FA_all,HR_all,'Color',colIdx(j,:))
            scatter(FA(i), HR(i),'g','Filled')
            title('ROC');
            xlabel('FA')
            ylabel('HR')
            xlim([0 1])
            ylim([0 1])
            
            pause(.1)
        end
        delete(h1)
    end
    hold off
end
pause;
clear all


%% B.3/6 - Optimal decision task
nSamples = 1000;
lambdaS = 10;
lambdaN = 5;

% Priors
pS = .2;

% Costs
lFA = 4;
lM  = 8;
lH  = 1;
lCR = 1;

[HR FA] = calcHR_FA(lambdaS,lambdaN,pS,nSamples,lFA,lM,lH,lCR,1);
fprintf('B.3\nSimulation\n\tHR: %.02f; FA: %.02f\n',HR,FA);
[HR FA] = calcHR_FA(lambdaS,lambdaN,pS,nSamples,lFA,lM,lH,lCR,0);
fprintf('Marginalization\n\tHR: %.02f; FA: %.02f\n\n',HR,FA);
pause;


%% B.4/6 - Altering pSource
pS = .5;
[HR FA] = calcHR_FA(lambdaS,lambdaN,pS,nSamples,lFA,lM,lH,lCR,1);
fprintf('B.4\nSimulation\n\tHR: %.02f; FA: %.02f\n',HR,FA);
[HR FA] = calcHR_FA(lambdaS,lambdaN,pS,nSamples,lFA,lM,lH,lCR,0);
fprintf('Marginalization\n\tHR: %.02f; FA: %.02f\n\n',HR,FA);
pause;


%% B.5/6 - Altering Costs
pS  = .2;
lFA = 8;
lM  = 4;
lH  = 1;
lCR = 1;
[HR FA] = calcHR_FA(lambdaS,lambdaN,pS,nSamples,lFA,lM,lH,lCR,1);
fprintf('B.4\nSimulation\n\tHR: %.02f; FA: %.02f\n',HR,FA);
[HR FA] = calcHR_FA(lambdaS,lambdaN,pS,nSamples,lFA,lM,lH,lCR,0);
fprintf('Marginalization\n\tHR: %.02f; FA: %.02f\n\n',HR,FA);











%%
function [iMAP, LS] = calcMAP_LS(x,y,MUx,SIGx)

% Get prior PDF
px = normpdf(x,MUx,SIGx);
for i = 1:length(y)
    % Likelihood
    py_x = normpdf(y(i),x,.05 + .2*x);
    
    % Normalization
    py = trapz(x,py_x .* px);
    
    % Posterior
    px_y = (py_x .* px) ./ py;
    
    [MAP(i), iMAP(i)] = findpeaks(px_y,'SortStr','descend','NPeaks',1);
    LS(i) = trapz(x,x .* px_y);
end


% Manually correct faulty peak finding at edge
iMAP(1) = 1;
iMAP(2) = 1;





%%
function [HR FA] = calcHR_FA(lambdaS,lambdaN,pS,nSamples,lFA,lM,lH,lCR,discrete)

pN = 1 - pS;

if discrete
    % Generate some observations given prior and likelihood
    cS = poissrnd(lambdaS, 1, nSamples*pS);
    cN = poissrnd(lambdaN, 1, nSamples*pN);
    c = [cS cN];
else
    c = 1:30;
end

% Index where there was actually signal
sIdx = zeros(1,nSamples);
sIdx(1:nSamples*pS) = 1;

% Estimate an internal measure given each observation - p(m|c):
pm_cS = poisspdf(c,lambdaS);
pm_cN = poisspdf(c,lambdaN);

% Estimated loss for each decision ("present" / "absent"):
% Formula = sum(p(m|c) * p(c) * L(c,c(m))) per decision
lossPresent = pm_cS .* pS .* lH + pm_cN .* pN .* lFA;
lossAbsent  = pm_cS .* pS .* lM + pm_cN .* pN .* lCR;

% A detection occurs when the loss from the signal is less than the
% loss from the noise
detect = lossPresent < lossAbsent;

if discrete
    % Hit rate = signal present, and you say signal is present
    HR = sum(sIdx & detect) / (nSamples * pS);

    % False alarm rate = signal absent, you say signal is present
    % (ie. when the loss
    FA = sum(~sIdx & detect) / (nSamples * pN);
else
    HR = trapz(c,detect.*pm_cS);
    FA = trapz(c,detect.*pm_cN);
end

    
    
    
    
    
    

    












