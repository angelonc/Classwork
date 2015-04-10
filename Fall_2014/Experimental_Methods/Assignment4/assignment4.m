function fitEstimates = assignment4
% poissonFOSFitExample
%
% This program demonstrates how to use fmincon to find
% the criterion for detection in a Poisson frequency of
% seeing model, when the fraction absorbed and dark noise
% are known
%
% 10/18/11  dhb  Wrote it.
% 10/18/14 CFA modified for multiple parameter estimates.


%% Clear and close
clear; close all;

%% Define parameters
fractionAbsorbed = .3;
darkNoise = .1;
detectCriterionAbsorbed = 7;
nTrials = 10;
testQuanta = logspace(1,3,20)';
theoryTestQuanta = logspace(1,3,100)';



%% Synthesize data based on the parameters
% First simulate performance
nTests = length(testQuanta);
absorbedQuanta = fractionAbsorbed*testQuanta;
nSee = zeros(size(testQuanta));
for i = 1:nTests
        actualAbsorbed = poissrnd(absorbedQuanta(i)+darkNoise,nTrials,1);
        nSee(i) = length(find(actualAbsorbed > detectCriterionAbsorbed));
end
percentSee = 100*nSee/nTrials;



%% Plot frequency of seeing curve against actual test intensity
% This depends on knowing the fraction absorbed.
figure; clf; hold on
set(gca,'FontName','Helvetica','FontSize',14);
plot(log10(testQuanta),percentSee,'ro','MarkerFaceColor','r','MarkerSize',10);
xlabel('Log10 Test Intensity','FontName','Helvetica','FontSize',18);
ylabel('Percent Seen','FontName','Helvetica','FontSize',18);
title(sprintf('Dark noise, %g; Detect criterion %g; Fraction absorbed %g',darkNoise,detectCriterionAbsorbed,fractionAbsorbed), ...
    'FontName','Helvetica','FontSize',14);



%% Put theoretical FOS curve through the data.  This is done knowing the parameters
% used to produce the simulated data, so it should be right up to effects of simulated
% measurement variability. 
theoryPercentSee = 100*poissonFractionSeen(theoryTestQuanta,detectCriterionAbsorbed,fractionAbsorbed,darkNoise);
plot(log10(theoryTestQuanta),theoryPercentSee,'r');



%% Now fit the FOS curve.  We use Matlab fmincon.  In this example, only
% the criterion is searched for, the variables fractionAbsorbed and darkNoise are
% assumed to be known and are used in the theoretical calculation.

% Standard options set up.  Sometimes one search algorithm works better than another.
% See the fmincon help for some information on what algorithms are available.
%
% Set 'Display' -> 'iter-detailed' to see more printout of what is happening.
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');

% Our parameter is the criterion absorbed.  That is, here we treat the fraction absorbed as known, and also
% assume that the dark noise is zero.
vlb = [0; 0; 0];
vub = [300; 10; .5];
params0 = [20; 1; .2];
fitEstimates = fmincon(@(params)poissonFitErrorFunction(params,testQuanta,percentSee,nTrials),params0,[],[],[],[],vlb,vub,[],options);
fprintf('Criterion estimate: %g (true value: %g)\n',fitEstimates(1),detectCriterionAbsorbed);
fprintf('Fraction absorption estimate: %g (true value: %g)\n',fitEstimates(2),fractionAbsorbed);
fprintf('Dark noise estimate: %g (true value: %g)\n', fitEstimates(3),darkNoise);

%% Add the fit curve to the plot in green.
fitPercentSee = 100*poissonFractionSeen(theoryTestQuanta,fitEstimates(1),fitEstimates(2),fitEstimates(3));
plot(log10(theoryTestQuanta),fitPercentSee,'g');
fileName = sprintf('poissonFitFOS_%g_%g_%g',darkNoise,detectCriterionAbsorbed,fractionAbsorbed);
% if (exist('savefig','file'))
%     savefig(fileName,gcf','pdf');
% else
%     saveas(gcf,fileName,'pdf');
% end


end

% f = poissonFitErrorFunction(criterion,testQuanta,fractionAbsorbed,darkNoise,percentSeen,nTrials)
% 
% The error function is the negative log likelihood of the data.
function f = poissonFitErrorFunction(params,testQuanta,percentSeen,nTrials)

criterion = params(1);
fractionAbsorbed = params(2);
darkNoise = params(3);
% Compute predicted fraction seen for each data point
predictedFractionSeen = poissonFractionSeen(testQuanta,criterion,fractionAbsorbed,darkNoise);
predictedFractionSeen(predictedFractionSeen == 0) = 0.0001;
predictedFractionSeen(predictedFractionSeen == 1) = 0.9999;

% Accumulate log likelihood for the entire data set by adding up the
% log likelihoods for each trial.
logLikelihood = 0;
fractionSeen = percentSeen/100;
numberSeen = round(nTrials*fractionSeen);
numberNotSeen = nTrials-numberSeen;
for i = 1:length(testQuanta)
    for j = 1:numberSeen(i)
        logLikelihood = logLikelihood + log10(predictedFractionSeen(i));
    end
    for j = 1:numberNotSeen(i)
        logLikelihood = logLikelihood + log10(1-predictedFractionSeen(i));
    end
end

% Since high likelihood is good, we take the negative to get the error function
% return value.
f = -logLikelihood;

end

% fractionSeen = poissonFractionSeen(testQuanta,criterion,fractionAbsorbed,darkNoise)
%
% Compute fraction seen under the poisson/criterion model.
%
% Inputs:
%   testQuanta: Vector, averagenumber of test quanta per flash.
%   criterion: Say seen if more than this number are absorbed.
%   fractionAbsorbed: fraction of test quanta that are absorbed and cause an isomerization.
%   darkNoise: number of spontaneous isomerizations per flash.  By this we mean teh
%     number of spontaneous isomerizations that would occur in the dark, in the area of
%     the retina corresponding to a flash over a duration corresponding to that of a flash.
%
% Outputs:
%   fractionSeen: Vector, fraction seen for test intensity specified in vector testQuanta.
%
% Because the predictions of the model are for discrete integer values of the
% criterion (can never absorb fractional photons), need to make them continuous
% by implementing a mixture of two nearest integer criteria, when passed criterion
% is not integer.

function fractionSeen = poissonFractionSeen(testQuanta,criterion,fractionAbsorbed,darkNoise)

criterionFloor = floor(criterion);
criterionCeil = ceil(criterion);
if (criterionFloor == criterionCeil)
    fractionSeen = (1-poisscdf(criterionCeil,fractionAbsorbed*testQuanta+darkNoise));
else
    fractionSeen = (criterion-criterionFloor)*(1-poisscdf(criterionCeil,fractionAbsorbed*testQuanta+darkNoise)) + ...
        (criterionCeil-criterion)*(1-poisscdf(criterionFloor,fractionAbsorbed*testQuanta+darkNoise));
end
end


