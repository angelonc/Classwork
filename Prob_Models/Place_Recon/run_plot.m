function run_plot(pv)
% Run reconstruction
if ~exist('Recon_Data.mat','file'); 
    [loc, recon, post, pv_recon] = place_recon(pv);
    save('Recon_Data.mat', 'loc','recon','post','pv_recon');
else
    load Recon_Data.mat
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
    imagesc(posterior);
    
    [~,argmax_idx] = max(post{i}(:));
    [x y] = ind2sub(size(post{i}),argmax_idx);
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
    hold off
    pause(dt/10);
end
figure(2);
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
        
        plot(pv_recon(1:i,2),pv_recon(1:i,1),'-b');
        scatter(pv_recon(i,2),pv_recon(i,1),'b','filled');
        
    else
        plot(loc(i-tail:i,2),loc(i-tail:i,1),'-g');
        scatter(loc(i,2),loc(i,1),'g','filled');
        
        plot(pv_recon(i-tail:i,2),pv_recon(i-tail:i,1),'-b');
        scatter(pv_recon(i,2),pv_recon(i,1),'b','filled');
    end
    hold off
end






