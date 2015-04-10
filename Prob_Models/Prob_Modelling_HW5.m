%% Probabilistic Modeling HW5
function Prob_Modelling_HW5
close all

%% A.1 - Bayesian Estimation

% Likelihood params:
mut = 70;
sigt = 2.5;

% Prior params:
mup = 50;
sigp = 5;

n = 1000;

% Set x range, simulate some neural measures
x = linspace(1,100,n);
m = normrnd(mut,sigt,n,1);

% Prior
px = normpdf(x,mup,sigp);

% Inference
for i = 1:n
    % Likelihood
    pm_x = normpdf(m(i),x,sigt);
    % Posterior
    px_m = (pm_x .* px) / trapz(x,pm_x .* px);
    % BLS estimate
    x_hat(i) = trapz(x,px_m .* x);
end

% Normalize and plot
plotx = 1:100;
figure(1)
hold on
set(gca,'FontSize',16);
[N,X] = hist(x_hat,plotx);
bar(X,N/trapz(X,N))
title('Predicted and Analytical p(xhat|x)');
xlim([50 80])
xlabel('Line Length (x in mm)');
ylabel('Probability');
pause;


% Derive p(x_hat|x) - see stocker & simoncelli, 2006
% Essentially it is the posterior in that condition?:
pm_x_ideal = normpdf(plotx,mut,sigt);
px_ideal = normpdf(plotx,mup,sigp);
pxh_x = (px_ideal .* pm_x_ideal) ./ trapz(plotx,px_ideal .* pm_x_ideal);
plot(plotx, pxh_x, 'r', 'LineWidth', 3);
hold off

%% A.2

sigr = 1.5;

% Simulate data
nTrials = 100;
xr = 55:2.5:85;

x = 0:.5:100;
px = normpdf(x,mup,sigp);

for i = 1:length(xr)
    % Simulate measures for the test and reference stimuli
    m_test = normrnd(mut,sigt,nTrials,1);
    m_ref = normrnd(xr(i),sigr,nTrials,1);
    
    % Now generate estimates for each trial based on the measures
    for j = 1:nTrials
        % Test estimate
        pmt_t = normpdf(m_test(j),x,sigt);
        pt_mt = (pmt_t .* px) ./ trapz(x,pmt_t .* px);
        xh_t = trapz(x,pt_mt .* x);
        
        % Ref estimate
        pmr_r = normpdf(m_ref(j),x,sigr);
        pr_mr = (pmr_r .* px) ./ trapz(x,pmr_r .* px);
        xh_r = trapz(x,pr_mr .* x);
        
        resp(i,j) = xh_r > xh_t;
    end
end     

% Solve for p(ref>test)
p_refgtest = bayesEstimator(mup,sigp,mut,sigt,sigr,x,xr);

% Plot psychometrics
figure(2)
set(gca,'FontSize',16);
hold on
h2(2) = plot(xr,p_refgtest,'r','LineWidth',3);
h2(1) = errorbar(xr,mean(resp,2),std(resp,0,2)/sqrt(nTrials),'b.', ...
         'MarkerSize',30);
xlim([55 85]);
title('2AFC Psychometric Curve');
xlabel('Reference Value (mm)');
ylabel('p(Ref.>Test)');
% Plot PSE
h2(3) = plot([xr(6) xr(6)],[0 p_refgtest(6)], 'r--');
legend(h2, 'Empirical', 'Actual', 'PSE','Location','SE');
pause;


%% B - MODEL FITTING
% Data varied the contrast of the text 
clear all
path = 'Misc/ex5Data.mat';
if ~exist(path,'file');
    path = input('Type data path: ','s');
end
load(path);

%% B.1 - Gaussian CDF Fit

% Parameters
dat = [D1(:,2) D2(:,2) D3(:,2)];
options = optimset('fminunc');
options = optimset(options,'Diagnostics','off','Display','off', ...
                   'LargeScale','off','Algorithm','active-set');
func = 'normCDF';
parStr = {'mu', 'sigma'};
params0 = [55; 10];

% Plot psychometric curves and fit curves
figure(3)
colIdx = [1 0 0; 0 1 0; 0 0 1];
fprintf('PARAMETER ETIMATES FOR NORMCDF MODEL:\n');
for i = 1:size(dat,2)
    % Plot empirical
    subplot(1,3,i)
    set(gca,'FontSize',16);
    hold on
    h3(1) = scatter(xr,dat(:,i)/n,30,colIdx(i,:));
    xlim([45 85]);
    title(sprintf('Contrast %g',i));
    xlabel('Reference Length (mm)');
    ylabel('p(Ref.>Test)');
    
    % Find fit
    [parFit(i,:)] = fminunc(@(params)fitLLFunc(func,params,dat(:,i),xr,n), ...
                            params0,options);

    cdf_nll(i,:) = fitLLFunc(func,parFit(i,:),dat(:,i),xr,n);
    PSE(i) = norminv(.5,parFit(i,1),parFit(i,2));
    dThresh(i) = norminv(.75,parFit(i,1),parFit(i,2)) - PSE(i);
    
    % Plot fit and PSE and thresholds
    h3(2) = plot(xr,normcdf(xr,parFit(i,1),parFit(i,2)),...
         'Color', colIdx(i,:),'LineWidth',3);
    h3(5) = plot([PSE(i) PSE(i)], [0 .5],'--','Color',colIdx(i,:));
    
    fprintf('\tC%g: %s = %.2f, %s = %.2f --- PSE = %.2f, D_thresh = %.2f\n', ...
            i,parStr{1},parFit(i,1),parStr{2},parFit(i,2), ...
            PSE(i),dThresh(i));
    
    hold off
end

pause;



%% B.2 - Standard SDT Fit
% Fit data simultaneously
x = 0:.5:100;
params0 = [55 55 2 2 2];
parFit = fminunc(@(params)fitLLFunc_SDT(func,params,dat,x,xr,n), ...
                            params0,options);
sdt_nll = fitLLFunc_SDT(func,parFit,dat,x,xr,n);

pars = [70 parFit(3); ...
        parFit(1) parFit(4); ...
        parFit(2) parFit(5)];

% Plot on previous figure
figure(3)
fprintf('\nPARAMETER ESTIMATES FOR SDT MODEL:\n');
colIdx = [.5 0 0; 0 .5 0; 0 0 .5];
for i = 1:size(dat,2)
    subplot(1,3,i)
    hold on
    set(gca,'FontSize',16);
    h3(3) = plot(xr,normcdf(xr,pars(i,1),pars(i,2)),'--', ...
                 'Color',colIdx(i,:),'LineWidth',3);
    fprintf('\tC%g: mu = %.2f, sigma = %.2f\n', i, pars(i,1), pars(i,2));
    %legend(h3,'Data','ML Fit','SDT Fit','PSE','Location','NW');
    hold off
end

pause;



%% B.3 - Bayes Model
% Fitting
params0 = [70 70 4 4 4 30 1];
parFit = fminunc(@(params)fitLLFunc_Bayes(func,params,dat,x,xr,n), ...
                 params0,options);
bayes_nll = fitLLFunc_Bayes(func,parFit,dat,x,xr,n);

pars = [70 parFit(3); ...
        parFit(1) parFit(4); ...
        parFit(2) parFit(5); ...
        parFit(6) parFit(7)];

% Plot on previous figure
figure(3)
colIdx = [1 1 0; 0 1 1; 1 0 1];
fprintf('\nPARAMETER ESTIMATES FOR BAYES MODEL:\n');
for i = 1:size(dat,2)
        subplot(1,3,i)
    hold on
    set(gca,'FontSize',16);
    h3(4) = plot(xr,bayesEstimator(pars(4,1), pars(4,2),pars(i,1), pars(i,2),pars(1,2),x,xr),':', ...
                 'Color',colIdx(i,:),'LineWidth',3);
    ylim([0 1]);
    fprintf('\tC%g: mu = %.2f, sigma = %.2f\n', i, pars(i,1), pars(i,2));
    legend(h3,'Data','ML Fit','SDT Fit','Bayes Fit','PSE','Location','NW');
    hold off
end
fprintf('\tPrior: mu = %.2f, sigma = %.2f\n', pars(4,1), pars(4, ...
                                                  2));

pause;

%% B.4 - Model Comparison
% Coin flip
likelihood = 0;
% For each condition
for c = 1:3
% For each reference value
    for i = 1:length(xr)
        % How many times did the data say 'ref>test'?
        K = dat(i,c);
        
        % Get likelihood for 'ref>test' (I think I can add here because log products become sums?)
        for j = 1:K
            likelihood = likelihood + log10(.5);
        end
        
        % Multiplied by likelihood for 'test>ref'
        for j = 1:n-K
            likelihood = likelihood + log10(1 - .5);
        end
    end
end

% Get negative
coin_nll = -1 * likelihood;




% Explaining itself
% Coin flip
likelihood = 0;
% For each condition
for c = 1:3
% For each reference value
    for i = 1:length(xr)
        % How many times did the data say 'ref>test'?
        K = dat(i,c);
        
        % Get likelihood for 'ref>test' (I think I can add here because log products become sums?)
        for j = 1:K
            likelihood = likelihood + log10(K/n);
        end
        
        % Multiplied by likelihood for 'test>ref'
        for j = 1:n-K
            likelihood = likelihood + log10(1 - (K/n));
        end
    end
end

% Get negative
emp_nll = -1 * likelihood;

normLL(1) = (bayes_nll - coin_nll) / (emp_nll - coin_nll)
normLL(2) = (sdt_nll - coin_nll) / (emp_nll - coin_nll)
normLL(3) = (sum(cdf_nll) - coin_nll) / (emp_nll - coin_nll)

figure(4)
set(gca,'FontSize',16);
b1 = bar(normLL);
set(gca,'XTickLabel',{'Bayes','SDT','CDF'});







































% Set up log-likelihood function to give to fmincon
function nll = fitLLFunc(func,params,D,xr,n)

% Get predicted probabilities for each reference value
% NormCDF model:
sigma = params(2);
mu = params(1);
modelEst = normcdf(xr,mu,sigma);


% Exception for 0 and 1
modelEst(modelEst == 0) = .001;
modelEst(modelEst == 1) = .999;

% Calculate log likelihoods (product) of each prediction
nll = nllFunc(size(modelEst,1), size(modelEst,2), n, modelEst, D);




% Set up log-likelihood function to give to fmincon
function nll = fitLLFunc_SDT(func,params,D,x,xr,n)

% Get predicted probabilities for each reference value
% SDT model:
mu = [70 params(1) params(2)];
sig = [params(3) params(4) params(5)];


% For each condition
for i = 1:length(mu)
    % Find cum area of each reference in the test distribution
    pm_t = normpdf(x,mu(i),sig(i));
    b = cumtrapz(x,pm_t);

    for j = 1:length(xr)
        % Reference
        pm_r = normpdf(x,xr(j),sig(1));
        est(i,j) = trapz(x,pm_r.*b);
    end
end




% Exception for 0 and 1
est(est==0) = .0001;
est(est>.9999) = .9999;



% Calculate log likelihoods (product) of each prediction
nll = nllFunc(size(est,1), size(est,2), n, est, D);





% Set up log-likelihood function to give to fmincon
function nll = fitLLFunc_Bayes(func,params,D,x,xr,n)

% Get predicted probabilities for each reference value
% Bayes model:
mu = [70 params(1) params(2)];
sig = [params(3) params(4) params(5)];
mup = params(6);
sigp = params(7);

px = normpdf(x,mup,sigp);
% For each condition
for i = 1:length(mu)
    % Solve for p(ref>test)
    est(i,:) = bayesEstimator(mup,sigp,mu(i),sig(i),sig(1),x,xr);
end

% Exception for 0 and 1
est(est==0) = .0001;
est(est>.9999) = .9999;

nll = nllFunc(size(est,1), size(est,2), n, est, D);




function est = bayesEstimator(mup,sigp,mut,sigt,sigr,x,xr)

% Solve for p(ref>test)

% Prior
px = normpdf(x,mup,sigp);

% Posterior of test
pxt_t = (normpdf(x,mut,sigt) .* px) ./ trapz(x,normpdf(x,mut,sigt) .* px);
b = cumtrapz(x,pxt_t);

px = normpdf(x,mup,sigp);
for i = 1:length(x)
    % Posterior of reference
    pxr_r = (normpdf(x,x(i),sigr) .* px) ./ trapz(x,normpdf(x,x(i),sigr) ...
                                                   .* px);
    p_refgtest(i) = trapz(x,pxr_r.*b);
end

% Get est for each xr value
for i = 1:length(xr)
    est(1,i) = p_refgtest(x==xr(i));
end


function nll = nllFunc(nConds, nVals, n, est, D)

% Calculate log likelihoods (product) of each prediction
likelihood = 0;
% For each condition
for c = 1:nConds
% For each reference value
    for i = 1:nVals
        % How many times did the data say 'ref>test'?
        K = D(i,c);
        
        % Get likelihood for 'ref>test' (I think I can add here because log products become sums?)
        for j = 1:K
            likelihood = likelihood + log10(est(c,i));
        end
        
        % Multiplied by likelihood for 'test>ref'
        for j = 1:n-K
            likelihood = likelihood + log10(1 - est(c,i));
        end
    end
end

% Get negative
nll = -1 * likelihood;









