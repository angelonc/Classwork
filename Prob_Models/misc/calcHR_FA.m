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