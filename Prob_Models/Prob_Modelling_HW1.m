%% PROBABALISTIC MODELING -- HW 1
function Prob_Modelling_HW1

close all


%% A1
fprintf('A.1) Bernoulli''s Urn\n');
fprintf('\n\tA.1.1) Red ball followed by 2 whites {''RWW''} for Urn A and B\n')
P_SWW_A = .6 * .4 * .4;
fprintf('\t\t**Prob of drawing RWW from urn A: %g\n', P_SWW_A);
P_SWW_B = .4 * .6 * .6;
fprintf('\t\t**Prob of drawing RWW from urn B: %g\n', P_SWW_B);
pause;

fprintf(['\n\tA.1.2) Drawing exactly two white balls in a draw of three: ' ...
         '{RWW, WRW, WWR}\n'])
P_WWX_A = P_SWW_A * 3;
fprintf('\t\t**P(exactly two whites from urn A): %g\n',P_WWX_A);
P_WWX_B = P_SWW_B * 3;
fprintf('\t\t**P(exactly two whites from urn B): %g\n', ...
        P_WWX_B);
pause;

fprintf(['\n\tA.1.3) Drawing at least two white balls in a draw of three: ' ...
         '{RWW, WRW, WWR, WWW}\n'])
P_WWW_A = .4^3;
P_WWW_B = .6^3;
fprintf('\t\t**P(at least two whites from urn A): %g\n',...
        P_WWX_A + P_WWW_A);
fprintf('\t\t**P(at least two whites from urn B): %g\n', ...
        P_WWX_B + P_WWW_B);
pause;

fprintf(['\n\tA.1.4) Say three balls were drawn from an urn, one red, ' ...
         'two white (Y = RWW).\n\tWhat is P(A|Y), P(B|Y)?\n'])
P_U = .5;
P_Y = (P_WWX_A * P_U) + (P_WWX_B * P_U);
P_YA = P_WWX_A;
P_YB = P_WWX_B;
fprintf('\t\t**P(A|Y): %g\n', (P_YA*P_U)/P_Y);
fprintf('\t\t**P(B|Y): %g\n', (P_YB*P_U)/P_Y);
pause;



%% A2
fprintf('\n\nA.2) Sampling without replacement.\n');
fprintf(['w = # white balls, r = # red balls. Two balls drawn without ' ...
         'replacement. Prove P(r1) = P(r2)\n']);
fprintf('\tt = r + w\n');
fprintf(['\tPossible draws = [rr rw wr ww]\n\tOnly care about ' ...
         'occurrances of r in either position:\n\tP(rr) = r/(t) * ' ...
         '(r-1)/(t-1) = (r^2-r)/(t^2+t)\n\tP(rw) = r/(t) * w/(t-1) = rw/(t^2+t)\n\tP(wr) = ' ...
         'w/(t) * r/(t-1) = rw/(t^2+t)\n\n']);
fprintf(['\t**P(r1) = P(rr) + P(rw) = (r^2 + rw - r)/(t^2 - t)\n\t**P(r2) ' ...
         '= P(rr) + P(wr) = (r^2 + rw - r)/(t^2 - t)\n']);
fprintf('\t**P(r1) = P(r2)\n');
pause



%% A3
fprintf('\n\nA.3) Conditional Probabilities - Fred, Alf, Bob\n')
fprintf('State space = [BAF ABF BFA AFB FBA FAB]\n');
fprintf(['\n\tA.3.1) P(Fred older than Bob)?\n\t\tTrue states = ' ...
         'BAF, ABF, BFA\n\t\t**P(F>B) = (# states F > B / total states) ' ...
         '= 3/6 = .5\n']);
fprintf(['\tA.3.2) P(Fred older than Bob and Alf)?\n\t\tTrue states = ' ...
         'BAF, ABF\n\t\t**P(F>B,A) = 2/6 = .33\n']);
fprintf(['\tA.3.3) P(Fred older than Bob | Fred older than Alf)?\n\t\tNew state space = ' ...
         'BAF, ABF, AFB\n\t\tTrue states = BAF, FAB\n\t\t**P(F>B|F>A) = ' ...
         '2/3 = .66\n']);
pause



%% A4
fprintf('\n\nA.4) Monty Hall Scenario\n');
fprintf(['\n\tA.4.1) You choose door 1, host chooses door 3, should ' ...
         'you switch choice to door 2?\n']);
fprintf(['\tLets set up some priors - in different prize conds,' ...
         ' what is P host will choose door 3?\n\t\tP(H3|D1) ' ...
         '= .5 (can''t choose D1, so .5 between D2 and D3)\n\t\tP(H3|D2) ' ...
         '= 1 (has to choose D3, b/c you chose D1 and the ' ...
         'prize is in D2)\n\t\tP(H3|D3) = 0 (can''t choose D3 if ' ...
         'prize is there)\n']);
priors = [.5 1 0];
fprintf(['\tWhat is prob of host choosing 3 given I chose 1?\n\t\tP(H3) ' ...
         '= .5\n\tWhat is prob of prize being in each door?\n\t\tP(D1) = P(D2) = P(D3) = .33\n'])
P_H3 = .5;
P_D = 1/3;
likelihood = P_D / P_H3;
fprintf(['\tNow, find prize probabilities given host choice of door ' ...
         '3:\n']);
post = priors * likelihood;
fprintf('\t\tP(D1|H3) = (P(H3|D1) * P(D1)) / P(H3) = %g\n', post(1));
fprintf('\t\tP(D2|H3) = (P(H3|D2) * P(D2)) / P(H3) = %g\n', post(2));
fprintf('\t\tP(D3|H3) = (P(H3|D3) * P(D3)) / P(H3) = %g\n', post(3));
[m,i] = max(post);
fprintf(['\t**Switch to door %g, more likely the prize is ' ...
         'there.\n'], i);
fprintf('\n\tA.4.2) Earthquake opened door 3.\n')
fprintf(['\t**Doesn''t matter if you switch, the earthquake ' ...
         'has no knowledge, priors = 1/3.\n']);
pause


%% A5
fprintf('\n\nA.5) Probability transformation\n');

% A.5.1
fprintf('\n\tA.5.1\n');
fprintf(['\n\tp''x(x-hat) = 0 = p''y(y-hat)\n\tx-hat = g(y-hat), ' ...
         'so y-hat = g^-1(x-hat)\n\tpy(y-hat)=px(g(y))*|g''(y-hat)|' ...
         '\n']);
fprintf(['\n\tp''y(y-hat) = (px(g(y-hat)))'' * (|g''(y-hat)|)'' ' ...
         '+ px(g(y-hat))*(|g''(y-hat)|)''\n']);
fprintf(['\n\tp''y(y-hat) = p''x(g(y-hat)) * (|g''(y-hat)|)'' * ' ...
         'g''(y-hat) + px(g(y-hat)) * (|g''(y-hat)|)''\n']);
fprintf('\n\tGiven y-hat = g^-1(x-hat)...');
fprintf(['\n\tp''y(y-hat) = p''x(x-hat) * g''(g^-1(x-hat)) ' ...
         '* (|g''(g^-1(x-hat))|)'' + px(x-hat) * (|g''(g^-1(x-hat))|' ...
         ')''']);
fprintf(['\n\tp''y(y-hat) = 0 * g''(g^-1(x-hat)) * (|g''(g^-1(x-hat))|)'' ' ...
         '+ px(x-hat) * (|g''(g^-1(x-hat))|)''']);
fprintf(['\n\n\tp''y(y-hat) = px(x-hat) * g''(g^-1(x-hat)) * g"(g^-1(x-hat)) ' ...
         '/ |g''(g^-1(x-hat))| = 0']);
fprintf('\n\t**Therefore max(px(x)) != max(py(y))');

fprintf('\n\n\tSo, the extrema of px will not be equal to the extrema of py.\n');
pause

% A.5.2
% Some parameters and variables
warning off
mu = 5;
sig = 1;
y_range = [-2 2];
x_range = [-6 10];
y = linspace(y_range(1),y_range(2),1000);
x = linspace(x_range(1),x_range(2),1000);

% Transformation function and its derivative (derivative of e^x = e^x)
xy = exp(y+2);
yx = log(x)-2;
dxy = xy;

% Probability density of x
px = normpdf(x,mu,sig);

% Find pdf for y given transformation function
py = normpdf(xy,mu,sig) .* abs(dxy);

% Get peaks of each distribution
[pkx,ppx] = findpeaks(px);
[pky,ppy] = findpeaks(py);

% Find peak of p(x) in y and peak of p(y) in x
pfx = normpdf(pky,mu,sig);
pfy = exp(y(ppx)+2);



% Plotting
fprintf('\n\tPlotting A.5.2...\n');
subplot(1,4,1);
set(gca, 'FontSize', 20);
hold on
plot(x,px,'r','LineWidth',3);
plot([x(ppx) x(ppx)], [0 pkx], 'r');
scatter(x(ppx),pkx,'r','filled','MarkerEdgeColor','k');
plot([xy(ppy) xy(ppy)], [0 .39], 'b--');
scatter(xy(ppy),.39,'b');
hold on
xlim(x_range)
ylim([0 2.25])
title('PDF p(x)');
xlabel('x');
ylabel('Probability Density x');

subplot(1,4,2)
set(gca, 'FontSize', 20);
hold on
plot(y,xy,'g','LineWidth',3);
scatter(y(ppy),xy(ppy),'b');
plot([y(ppy) y(ppy)], [-8 xy(ppy)], 'b--');
hold off
title('Transformation Function x(y)')
xlim(y_range)
ylim(x_range)
xlabel('y');
ylabel('x(y)');

subplot(1,4,3)
set(gca, 'FontSize', 20);
hold on
plot(x,yx,'g','LineWidth',3);
plot([x(ppx) x(ppx)], [-8 yx(ppx)], 'r--');
scatter(x(ppx),yx(ppx),'r');
hold off
title('Transformation Function y(x)')
xlim(x_range)
ylim(y_range)
xlabel('x');
ylabel('y(x)');

subplot(1,4,4)
set(gca, 'FontSize', 20);
hold on
plot(y,py,'b','LineWidth',3)
plot([y(ppy) y(ppy)], [0 pky], 'b')
scatter(y(ppy),pky,'b','filled','MarkerEdgeColor','k');
plot([yx(ppx) yx(ppx)], [0 1.98], 'r--');
scatter(yx(ppx),1.98,'r');
hold on
xlim(y_range)
ylim([0 2.25])
title('PDF p(y)');
xlabel('y');
ylabel('Probability Density y');
pause

clear all


%% A6
fprintf('\n\nA.6) - Inference\n\n');
x_range = [0 30];
y_range = [0 20];

% Simulate
x_sim = linspace(x_range(1),x_range(2),1000);
for i = 1:length(x_sim)
    y_sim(i) = normrnd(x_sim(i), .5 + (.2 * x_sim(i)));
end
figure
set(gca, 'FontSize', 20);
scatter(x_sim,y_sim,5,'filled');
title('A.6 Simulation');
xlabel('Part Length x');
ylabel('Measurement y');
fprintf('\tSimulated data set given the parameters\n');


% A.6.1
y = linspace(y_range(1),y_range(2),1000);
x = linspace(x_range(1),x_range(2),1000);
colorIdx = [1 0 0; 1 .5 0];
xi = [5 10];
pause

fprintf('\n\tPlotting A.6.1...\n');
figure
set(gca, 'FontSize', 20);
hold on
for i = 1:2
    % Conditional probability function:
    Sigma = .05 + .2 * (xi(i));
    Mu = xi(i);
    py_x5_10 = 1./(Sigma.*sqrt(2.*pi)) ...
        .* exp(-1.*((y-Mu).^2)./(2.*Sigma.^2));
    % py_x5_10 = normpdf(y, xi(i), .05 + (.2.*xi(i)));
    ph(i) = plot(y, py_x5_10, 'Color', colorIdx(i,:), 'LineWidth', 3);
    areas(i) = trapz(y,py_x5_10);
    fprintf('\tArea of PDF(Y|X=%g) = %g\n', xi(i), areas(i));
end
hold off
title('Forward PDF of Y=y|X = 5,10mm');
xlabel('Measurement y (mm)');
ylabel('Probability Density');
legend(ph, 'Xi = 5', 'Xi = 10');

pause

% A.6.2
fprintf('\n\tPlotting A.6.2...\n')
yi = 7.5;
Sigma = .05 + .2 .* x;
Mu = x;
py75_x = 1./(Sigma.*sqrt(2.*pi)) ...
        .* exp(-1.*((yi-Mu).^2)./(2.*Sigma.^2));
% py75_x = normpdf(yi, Mu, Sigma);

figure(100)
hold on
set(gca,'FontSize',20);
h(1) = plot(x,py75_x, 'b', 'LineWidth', 3);
title('Likelihood, Prior, Posterior PDFs');
xlim([0 x_range(2)]);
xlabel('Part Length x (mm)');
ylabel('Probability Density');

% A.6.3
fprintf('\n\tPlotting A.6.3...\n')
[ML, MLi] = findpeaks(py75_x);
plot([x(MLi) x(MLi)], [0 ML], 'b');
text(x(MLi)+.005,ML+.005,sprintf('ML = %.02f', x(MLi)));
pause

% A.6.4
fprintf('\n\tA.6.4\n');
Mu = 10;
Sigma = sqrt(8);

fprintf('\n\tSimulation for A.6.4 (no bounds)\n');
figure
set(gca,'FontSize',20);
its = 5000;
% Simulate
for i = 1:its
    x_sim2(i) = normrnd(Mu,Sigma);
    y_sim2(i) = normrnd(x_sim2(i),.05 + .2*x_sim2(i));
end
scatter(x_sim2,y_sim2,10,'filled');
title('6.4 Simulation');
xlabel('Part Length x');
ylabel('Measurement y');
pause

% Prior
px = normpdf(x,Mu,Sigma);

% Posterior
posterior = (py75_x .* px) ./ (trapz(x,py75_x .* px));
[MAP, iMAP] = findpeaks(posterior);

fprintf('\tPlotting A.6.4...\n');
figure(100)
h(2) = plot(x,px,'r','LineWidth',3);
h(3) = plot(x,posterior,'Color',[1 0 1],'LineWidth',3);
plot([x(iMAP) x(iMAP)], [0 MAP], 'Color', [1 0 1]);
text(x(iMAP)+.005,MAP+.005,sprintf('MAP = %.02f', x(iMAP)));
legend(h,'Conditional', 'Prior', 'Posterior')
hold off
pause

fprintf(['\n\n\tA.6.5: The ML estimate doesn''t take prior information ' ...
         'into account, while the MAP does.\n\tThis causes the MAP ' ...
         'estimate to be shifted closer to the prior p(x).\n']);







    
