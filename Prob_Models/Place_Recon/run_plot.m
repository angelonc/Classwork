function run_plot(pv,human)
% Run reconstruction

if ~human
    if ~exist('Recon_Data.mat','file'); 
        [loc, recon, post, recon_nocc, post_nocc, pv_recon] = place_recon(pv);
        save('Recon_Data.mat', 'loc','recon','post','pv_recon','recon_nocc');
    else
        load Recon_Data.mat
    end
else
    if ~exist('Recon_Data_Human.mat','file'); 
        [loc, recon, post, recon_nocc, post_nocc, pv_recon, l, u] = place_recon_hu(pv);
        save('Recon_Data_Human.mat','loc','recon','post', ...
             'pv_recon','recon_nocc','post_nocc','l','u');
    else
        load Recon_Data_Human.mat 
    end
    
    % Make new locations in bin space
    loc = (loc - l) / (u - l) * 70;
end




n = length(loc);
dt = 1; % Time in seconds
tail = 10;

% Smoothing kernel
sigma = (3 / 2*sqrt(2*log(2)))/.5;
G_kern = fspecial('gaussian',[7*ceil(sigma) 7*ceil(sigma)], ...
                  sigma);

map = [linspace(0,.7,64)' zeros(64,1) zeros(64,1)];
colormap(map);

% Make the video
figure(1);
title('Bayesian Reconstruction (no CC)');
set(gca,'Color',[0 0 0]);
for i = 1:n
    cla
    % Bayes
    xlim([1 70]);
    ylim([1 70]);
    hold on
    tmp = post_nocc{i};
    tmp(isnan(tmp)) = 0;
    posterior = imfilter(tmp,G_kern);
    imagesc(posterior./sum(posterior(:)));
    if i <= tail
        
        plot(loc(1:i,2),loc(1:i,1),'-g');
        scatter(loc(i,2),loc(i,1),'g','filled');
        
        plot(recon_nocc(1:i,2),recon_nocc(1:i,1),'-r');
        scatter(recon_nocc(i,2),recon_nocc(i,1),'r','filled');
        
    else
        plot(loc(i-tail:i,2),loc(i-tail:i,1),'-g');
        scatter(loc(i,2),loc(i,1),'g','filled');
        
        plot(recon_nocc(i-tail:i,2),recon_nocc(i-tail:i,1),'-r');
        scatter(recon_nocc(i,2),recon_nocc(i,1),'r','filled');
    end
    set(gca,'FontSize',16);
    hold off
    pause(dt/10);
end
keyboard
% Make the video
map = [linspace(0,.7,64)' zeros(64,1) zeros(64,1)];
colormap(map);
figure(1);
title('Bayesian Reconstruction (Two Step)');
set(gca,'Color',[0 0 0]);
for i = 1:n
    cla
    % Bayes
    xlim([1 70]);
    ylim([1 70]);
    hold on
    tmp = post{i};
    tmp(isnan(tmp)) = 0;
    posterior = imfilter(tmp,G_kern);
    imagesc(posterior./sum(posterior(:)));
    if i <= tail
        
        plot(loc(1:i,2),loc(1:i,1),'-g');
        scatter(loc(i,2),loc(i,1),'g','filled');
        
        plot(recon(1:i,2),recon(1:i,1),'-r');
        scatter(recon(i,2),recon(i,1),'r','filled');
        
    else
        plot(loc(i-tail:i,2),loc(i-tail:i,1),'-g');
        scatter(loc(i,2),loc(i,1),'g','filled');
        
        plot(recon(i-tail:i,2),recon(i-tail:i,1),'-r');
        scatter(recon(i,2),recon(i,1),'r','filled');
    end
    set(gca,'FontSize',16);
    hold off
    pause(dt/10);
end
keyboard

figure(1);
title('Population Vector Reconstruction');
set(gca,'Color',[0 0 0]);
for i = 1:length(pv_recon)
    cla
    % Population Vector
    xlim([1 70]);
    ylim([1 70]);
    hold on
    
    [~,argmax_idx] = max(post{i}(:));
    [x y] = ind2sub(size(post{i}),argmax_idx);
    if i <= tail
        
        plot(loc(1:i,2),loc(1:i,1),'-g');
        scatter(loc(i,2),loc(i,1),'g','filled');
        
        plot(pv_recon(1:i,2),pv_recon(1:i,1),'-c');
        scatter(pv_recon(i,2),pv_recon(i,1),'c','filled');
        
    else
        plot(loc(i-tail:i,2),loc(i-tail:i,1),'-g');
        scatter(loc(i,2),loc(i,1),'g','filled');
        
        plot(pv_recon(i-tail:i,2),pv_recon(i-tail:i,1),'-c');
        scatter(pv_recon(i,2),pv_recon(i,1),'c','filled');
    end
    set(gca,'FontSize',16);
    hold off
    pause(dt/10);
end

keyboard
figure(2)
bayes_err = sqrt((recon(:,1) - loc(:,1)).^2 + ...
                 (recon(:,2) - loc(:,2)).^2);
bayes_nocc_err = sqrt((recon_nocc(:,1) - loc(:,1)).^2 + ...
                 (recon_nocc(:,2) - loc(:,2)).^2);
povec_err = sqrt((pv_recon(:,1) - loc(:,1)).^2 + ...
                 (pv_recon(:,2) - loc(:,2)).^2);
errorbar([mean(bayes_err/2) mean(bayes_nocc_err/2) mean(povec_err/2)], ...
         [std(bayes_err/2)/sqrt(n) std(bayes_nocc_err/2)/sqrt(n) std(bayes_err/2)/sqrt(n)],...
         'LineWidth',3);
title('Human Error Comparison');
xlabel('Reconstruction Model');
ylabel('Average Error (cm)');
ylim([7 25]);
set(gca,'XTick',[1 2 3]);
set(gca,'XTickLabels',{'Bayes (2 Step)', 'Bayes (no CC)', 'Pop. Vector'});
set(gca,'FontSize',16);

keyboard









