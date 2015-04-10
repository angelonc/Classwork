%% PROBABILISTIC MODELING HW 4
function Prob_Modelling_HW4
close all

%% A.1 Psychic Abilities
% 1) H0 = random guessing
pD_H0 = (1/6)^5 * (5/6)^5;

% 2) H1 = 50/50 shot per throw
pD_H1 = (1/2)^5 * (1/2)^5;

% 3) Likelihood ratio
bf = pD_H1 / pD_H0;

% 5)
post = [.01 .99; .5 .5];
for i = 1:2
    odds(i) = (pD_H1/pD_H0) * (post(i,1)/post(i,2));
end


%% A.2 Maximum-Likelihood Fits
% Load data
path = input('Type data path: ','s');
if isempty(path) path = 'Misc/ex4Data.mat'; end
load(path);

% 1) Plot:
figure;
set(gca,'FontSize',16);
h(1) = scatter(xr,D(:,2)/n,'b','Filled');
title('p(Reference>Test) x Reference Value')
xlabel('Reference Value');
ylabel('p(Reference>Test)');
hold on

% 2) ML fit:
options = optimset('fminunc');
options = optimset(options,'Diagnostics','off','Display','off', ...
                   'LargeScale','off','Algorithm','active-set');

% Do fit for linear stepwise function and normCDF
for i = 1:2
    if i == 1
        func = 'normCDF';
        parStr = {'mu', 'sigma'};
        params0 = [55; 10];
        vlb = [0; 1];
        vub = [140; 20];
    else
        func = 'swLinear';
        parStr = {'a','b'};
        params0 = [56; 80];
        vlb = [55; 55];
        vub = [85; 85];
    end
    
    %parFit(i,:) = fmincon(@(params)fitLLFunc(func,params,D,xr,n),params0,[],[],[],[],vlb,vub,[],options);
    
    [parFit(i,:)] = fminunc(@(params)fitLLFunc(func,params,D,xr,n), ...
                            params0,options);
    
    ll(i,:) = -fitLLFunc(func,parFit(i,:),D,xr,n);
    aic(i) = aicbic(ll(i,:),length(parFit(i,:)));
    

    
    fprintf('Estimated parameters for %s model: %s = %g, %s = %g; AIC = %g\n', ...
            func,parStr{1},parFit(i,1),parStr{2},parFit(i,2),aic(i));
end

h(2) = plot(xr,normcdf(xr,parFit(1,1),parFit(1,2)),'r','LineWidth',2);
h(3) = plot(xr,linearSW(xr,parFit(2,1),parFit(2,2)),'g', ...
            'LineWidth',2);
legend(h,'Data','normCDF','SW Linear','Location','SE');
hold off












% Set up log-likelihood function to give to fmincon
function nll = fitLLFunc(func,params,D,xr,n)

% Get predicted probabilities for each reference value
if strcmp(func, 'normCDF');
    % NormCDF model:
    sigma = params(2);
    mu = params(1);
    modelEst = normcdf(xr,mu,sigma);
elseif strcmp(func, 'swLinear')
    % Linear step-wise model:
    a = params(1);
    b = params(2);
    modelEst = linearSW(xr,a,b);
end

% Exception for 0 and 1
modelEst(modelEst == 0) = .001;
modelEst(modelEst == 1) = .999;

% Calculate log likelihoods (product) of each prediction
likelihood = 1;
% For each reference value
for i = 1:length(xr)
    % How many times did the data say 'ref>test'?
    K = D(i,2);
    
    % Get likelihood for 'ref>test' (I think I can add here because log products become sums?)
    for j = 1:K
        likelihood = likelihood + log10(modelEst(i));
    end
    
    % Multiplied by likelihood for 'test>ref'
    for j = 1:n-K
        likelihood = likelihood + log10(1 - modelEst(i));
    end
end

% Get negative
nll = -1 * likelihood;



function y = linearSW(xr, a, b)

% Get slope and intercept
m = 1 / (b - a);
int = 0 - (m*a);

% Calculate discrete function
y = zeros(1,length(xr));
for i = 1:length(xr)
    if m > 0
        if xr(i) <= a
            y(i) = 0;
        elseif xr(i) >= b
            y(i) = 1;
        else
            y(i) = m*xr(i) + int;
        end
    elseif m < 0
        if xr(i) >= a
            y(i) = 0;
        elseif xr(i) <= b
            y(i) = 1;
        else
            y(i) = m*xr(i) + int;
        end
    else
        y(i) = 0;
    end
end


